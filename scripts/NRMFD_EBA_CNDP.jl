
using Pkg
Pkg.activate("CNDP")

using Printf
using LinearAlgebra
using Graphs, SparseArrays
using MutableNamedTuples
using Random
using Optim
using DelimitedFiles
using DataFrames, XLSX, CSV
using Distributions
using SimpleWeightedGraphs
using JuMP, Ipopt

include("../src/variables.jl")
include("../src/misc.jl")
include("../src/utils.jl")
include("../src/optimizer.jl")
include("../src/dataloader.jl")


function path_flow_counter(td, m_od)
    flow_path = Vector{Float64}()
    for r in 1:size(td.travel_demand)[1]
        for ss in 1:length(m_od[r])
            s = m_od[r][ss].D
            if td.travel_demand[r, s] == 0
                continue
            end
            for flow in m_od[r][ss].flows
                if flow > 0
                    push!(flow_path, flow)
                end
            end
        end
    end
    num_path = length(flow_path)
    return num_path, flow_path
end

function path_matrix(td, od_pairs, m_od, num_path)
    Delta = spzeros(td.number_of_links, num_path)
    Lambda = spzeros(length(od_pairs), num_path)
    timer = time()
    n_path_add = 0

    for r in 1:size(td.travel_demand)[1]
        for ss in 1:length(m_od[r])
            s = m_od[r][ss].D
            if td.travel_demand[r, s] == 0
                continue
            end
            @views tmp = m_od[r][ss]
            for (path, path_flow) in zip(tmp.paths, tmp.flows)
                if path_flow == 0
                    continue
                end
                n_path_add += 1
                for ind in path
                    Delta[ind, n_path_add] = 1
                end
                for j in findall(x -> x == [r, s], od_pairs)
                    Lambda[j, n_path_add] = 1
                end
            end
        end
    end
    return vcat(Delta, Lambda)
end

function path_split(td, flow, ep, Delta, Lambda, num_path, threshold)
    dtv = dtdv(td, flow, ep)
    dvy = Array{Float64}(undef, size(flow[flow.>0])[1], size(ep.y)[1])
    Delta = Delta[flow.>0, :]
    full_matrix = vcat(Delta, Lambda)
    old = num_path
    insetA = Vector{Int64}()
    while true 
        # if the threshold is too small, the matrix A may be singular and throw error.
        # if the threshold is too large, the result may be not accuracy
        # this while loop is used to find an threshold
        # global insetA, insetB, insetC 
        try
            full_matrix = vcat(Delta, Lambda)
            if rank(full_matrix) < old 
                insetA = QR_decomposition(full_matrix, threshold)
                while true
                    if length(insetA) == old
                        threshold *= 10
                    else
                        break
                    end
                    insetA = QR_decomposition(full_matrix, threshold)
                end
                old = size(insetA)[1]
                @views Delta1 = Delta[:, insetA]
                @views Lambda1 = Lambda[:, insetA]
            else
                insetA = collect(1:num_path)
                @views Delta1 = Delta
                @views Lambda1 = Lambda
            end
            dth = Delta1' * (dtv) * Delta1
            dty = Delta1' * dtdy(td, flow, ep)
            
            A = vcat(hcat(dth, -Lambda1'),
            hcat(Lambda1, zeros(size(Lambda1)[1], size(Lambda1)[1])))
            # @printf("size of A is %4d x %4d, rank of A is %4d\n", size(A)[1], size(A)[2], rank(A))
            b = vcat(-dty, zeros(size(Lambda1)[1], size(dty)[2]))
            dvy = Delta1 * (A\collect(b))[1:size(dty)[1], :]
            println("working threshold: ", threshold)
            break
        catch
            threshold *= 10
            println("try new threshold: ", threshold)
            if threshold > 1; break; end
        end
    end
    insetB = setdiff(1:num_path, insetA)
    insetC = Vector{Int64}()
    return [insetA, insetB, insetC], dvy
end

function pos_step_cal(f, df, alpha0)
    alpha = alpha0
    tmp1, tmp2 = f, df 
    tmp3 = f - df * alpha0 
    if sum(tmp3 .< 0) > 0
        ind1 = findall(x -> x<0, tmp3)
        if !isempty(ind1)
            tmp1 = tmp1[ind1]
            tmp2 = tmp2[ind1]
            alpha = minimum(tmp1 ./ tmp2)
        end
    end
    return alpha
end

function update_positively(y, dy, alpha, threshold=1e-12)
    y -= dy * alpha
    y[y.< threshold] .= 0
    return y
end

function update_idx!(insets, pi_path, flow_path, full_matrix)
    @views insetA, insetB, insetC = insets[1], insets[2], insets[3]
    # update set 
    # update_idx!(insets, pi_path1, flow_path1, full_matrix)
    # case 1: if pi_path1[insetC] .== 0, 
    # then remove the path from the NUE set
    # move it to ELI set or ELD set 
    if sum(pi_path[insetC] .== 0) > 0
        ind1 = findall(x -> x == 0, pi_path[insetC])
        tmp = insetC[ind1]
        insetC = setdiff(insetC, tmp)
        flow_path[tmp] .= 0.001
        for itmp in tmp 
            tmp2 = union(insetA, itmp)
            if rank(full_matrix[:, tmp2]) > rank(full_matrix[:, insetA])
                insetA = tmp2 
            else 
                insetB = union(insetA, itmp)
            end
        end
    end

    # # case 2 & 3: flow_path1[insetA] .== 0
    while sum(flow_path[insetA] .== 0) > 0
        # find the alternative set R' 
        # every collumn records the coefficients of i in A to represent j in B 
        threshold = 1e-10
        coes = (full_matrix[:, insetA]' * full_matrix[:, insetA]) \ collect(full_matrix[:, insetA]') * full_matrix[:, insetB]
        coes[abs.(coes) .< threshold] .= 0
        insetRA = findall(x -> x == 0, flow_path[insetA])
        # println("check size : ", size(coes), " == ", length(insetA), ", ", length(insetB))
        # case 2: if it has no alternative paths (P4) 
        # or suppose it has alternative paths set b \in R' but 
        # (i) the flow of all ELD paths is 0 
        # (ii) the coefficient of this zero-flow ELI path in any b \in R' is positive 
        # this ELI path and all its alternative paths should be removed to the NUE set
        
        # case 3: replace it with an eligible ELD paths in R'
        # if there is b \in R' s.t. flow > 0, then choose it 
        # otherwise, choose a b \in R' s.t. the coefficient of this zero-flow ELI path in b is negative
        case2flag = false 
        case3flag = false 
        # println("Ra size: ", length(insetRA))
        rA = insetRA[1]
        # for rA in insetRA
        r = insetA[rA]
        indbb = findall(x -> x != 0, coes[rA, :])
        Rp = insetB[indbb]
        # println(coes[rA, indbb])
        cond4 = length(indbb) == 0
        if cond4 
            insetA = setdiff(insetA, r)
            insetC = union(insetC, r)
            pi_path[r] = 0.001
            pi_path[Rp] .= 0.001
            println("case 2")
        else 
            cond51 = (sum(flow_path[Rp] .== 0) == length(indbb))
            cond52 = (sum(coes[rA, indbb] .< 0) == 0)
            if cond51 & cond52 
                insetA = setdiff(insetA, r)
                insetB = setdiff(insetB, Rp)
                insetC = union(insetC, r, Rp)
                pi_path[r] = 0.001
                pi_path[Rp] .= 0.001
                println("case 2")
            elseif !cond51 
                indbb = findall(x -> x != 0, flow_path[Rp])[1]
                indb = Rp[indbb]
                # @printf("case 3a, rank old A is %4d\n", rank(full_matrix[:, insetA]))
                insetA = setdiff(insetA, r)
                insetA = union(insetA, indb)
                insetB = setdiff(insetB, indb)
                insetB = union(insetB, r)
                # @printf("case 3a, new flow is %.3f rank new A is %4d\n", flow_path1[indb], rank(full_matrix[:, insetA]))
            elseif !cond52 
                # case 3
                indbb = (indbb[findall(x -> x < 0, coes[rA, indbb])])[1]
                indb = Rp[indbb]
                insetA = setdiff(insetA, r)
                insetA = union(insetA, indb)
                insetB = setdiff(insetB, indb)
                insetB = union(insetB, r)
                flow_path[indb] = 0.001
                println("case 3b, new flow is ", flow_path[indb])
            else 
                @assert true println("edge case!!!!!! cond4 is ", cond4, " cond51 is ", cond51, " cond52 is ", cond52)
            end
        end
    end
    return insets
end

function QR_decomposition(A::SparseMatrixCSC{Float64,Int64}, threshold::Float64)
    Q, R, P = qr(collect(A), Val(true))
    # 找出 R 的对角线上绝对值大于阈值的元素
    cols = abs.(diag(R)) .> threshold
    # 计算矩阵的秩
    rank = sum(cols)

    indices = P[1:rank]

    return indices
end

function solve_direction(g, y)
sigma = 2
ny = length(y)

prob = Model(Ipopt.Optimizer)
set_silent(prob)
@variable(prob, d[i = 1:ny])
@variable(prob, eta)
@constraint(prob, y + d >= 0)
@constraint(prob, g' * d - eta <= 0)
@NLobjective(prob, Min, eta + sigma * sum(d[i] * d[i] for i in 1:ny))

set_start_value.(d, g)
optimize!(prob)

return - value.(d)

end

function NRMFD_EBA(; name = "SiouxFalls", jam=1, maxIter = 50, per = 0.1, seed = 1, dmean = 2., dstd = 4., due_acc = 1e-2, verbose = 2)
    # read data 
    results_root = "../results/CNDP/"
    td, flow, ep = load_tntp(name, per, seed, dmean, dstd, jam)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end

    # Initialization
    T = 0
    t0 = time()
    all_or_nothing!(m_od, flow, td, graph, link_dic)
    nii = 0
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        nii += 1
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        if nii > 50; break; end
    end

    T += time() - t0
    df0 = DataFrame(iter = 0, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    ys = [deepcopy(ep.y)]

    t0 = time()
    od_pairs0 = Vector{Vector}()
    # for (i, j) in zip(td.init_node, td.term_node)
    #     if td.travel_demand[i, j] > 0
    #         push!(od_pairs0, [i, j])
    #     end
    # end

    od_pairs = Vector{Vector}()
    travel_demand = Vector{Float64}()

    for i in 1:size(td.travel_demand)[1]
        for j in 1:size(td.travel_demand)[2]
            if td.travel_demand[i, j] > 0
                push!(od_pairs, [i, j])
                push!(travel_demand, td.travel_demand[i, j])
            end
        end
    end

    num_links = td.number_of_links
    println(length(od_pairs), " ", length(od_pairs0))

    # optimization 
    T += time() - t0

    for nIter in 1:maxIter
        t0 = time()

        # find the direction 
        # collect the path and flow information 
        num_path, flow_path = path_flow_counter(td, m_od)

        full_matrix = path_matrix(td, od_pairs, m_od, num_path)
        println("size of full_matrix is ", size(full_matrix[flow.>0, :]))
        @views Delta0 = full_matrix[1:num_links, :]
        @views Lambda0 = full_matrix[num_links+1:end, :]
        
        # split the path set into three subsets and calculate the direction 
        threshold = 1e-12 # initial threshold
        insets, dvy = path_split(td, flow, ep, Delta0, Lambda0, num_path, threshold)
        g0 = get_gradients(dvy, td, flow, ep, 0)
        # g = min.(g0, -ep.y)
        g = solve_direction(g0, ep.y)
        # g = g0
        if verbose > 1
            @printf("mean value of abs(g) is %.2f\n", mean(abs.(g)))
        end

        # decide the step size 
        decent_flag = true 
        # alpha0 = 0.15
        ep0   = deepcopy(ep)
        flow0 = deepcopy(flow)
        flow_path0 = deepcopy(flow_path)
        pi_path0 = zeros(num_path)
        alphas = Vector{Float64}([0])
        Fold = upper_obj(td, flow0, ep0)
        Fs = Vector{Float64}([Fold])
        p  = 0.01
        alpha0 = p * Fold / (g0' * g0)
        if verbose > 1
            println("Suggested step size: ", p * Fold / (g0' * g))
        end
        @views insetA, insetB, insetC = insets[1], insets[2], insets[3]

        while decent_flag
            flow_path1 = deepcopy(flow_path0)
            pi_path1 = deepcopy(pi_path0)
            flow1 = deepcopy(flow0)
            ep1 = deepcopy(ep0)
            
            inset = union(insetA, insetC)
            ninset = length(inset)
            @views Delta13 = Delta0[flow0.>0, inset]
            @views Lambda13 = Lambda0[:, inset]
            dtv = dtdv(td, flow0, ep0)
            dth = Delta13' * (dtv) * Delta13
            dty = Delta13' * dtdy(td, flow0, ep0)

            A = vcat(hcat(dth, Matrix{Float64}(I, ninset, ninset), -Lambda13'),
            hcat(diagm(pi_path0[inset]), diagm(flow_path0[inset]), zeros(ninset, length(od_pairs))),
            hcat(Lambda13, zeros(length(od_pairs), ninset + length(od_pairs))))
            b = vcat(-dty, zeros(ninset + length(od_pairs), size(dty)[2]))
            # @printf("size of A is %4d x %4d, rank of A is %4d\n", size(A)[1], size(A)[2], rank(A))
            
            dv = A \ collect(b)          
            dvy = Delta13 * dv[1:ninset, :]
            dvpi = dv[ninset+1:2*ninset, :]
            df = dv[1:ninset, :] * g 
            dpi = dvpi * g 
            dy = Delta13 * df
        
            # check the states
            alpha1 = pos_step_cal(flow_path1[inset], df, alpha0)
            alpha2 = pos_step_cal(pi_path1[inset], dpi, alpha0)
            alpha3 = pos_step_cal(flow1[flow1.>0], dy, alpha0)
            alpha = minimum([alpha1, alpha2, alpha3])
            if verbose > 0
                @printf("alpha is %.2f, which is minimum of ( %.2f, %.2f, %.2f )\n", alpha, alpha1, alpha2, alpha3)
            end

            eps = 1e-12
            flow_path1[inset] .= update_positively(flow_path1[inset], df, alpha, eps)
            pi_path1[inset] .= update_positively(pi_path1[inset], dpi, alpha, eps)
            flow1[flow1.>0] .= update_positively(flow1[flow1.>0], dy, alpha, eps)

            # # update set 
            # insets = update_idx!(insets, pi_path1, flow_path1, full_matrix)
            # case 1: if pi_path1[insetC] .== 0, 
            # then remove the path from the NUE set
            # move it to ELI set or ELD set 
            if sum(pi_path1[insetC] .== 0) > 0
                ind1 = findall(x -> x == 0, pi_path1[insetC])
                tmp = insetC[ind1]
                insetC = setdiff(insetC, tmp)
                flow_path1[tmp] .= 0.001
                for itmp in tmp 
                    tmp2 = union(insetA, itmp)
                    if rank(full_matrix[:, tmp2]) > rank(full_matrix[:, insetA])
                        insetA = tmp2 
                    else 
                        insetB = union(insetA, itmp)
                    end
                end
            end

            # # case 2 & 3: flow_path1[insetA] .== 0
            while sum(flow_path1[insetA] .== 0) > 0
                # find the alternative set R' 
                # every collumn records the coefficients of i in A to represent j in B 
                threshold = 1e-10
                coes = (full_matrix[:, insetA]' * full_matrix[:, insetA]) \ collect(full_matrix[:, insetA]') * full_matrix[:, insetB]
                coes[abs.(coes) .< threshold] .= 0
                insetRA = findall(x -> x == 0, flow_path1[insetA])
                # println("check size : ", size(coes), " == ", length(insetA), ", ", length(insetB))
                # case 2: if it has no alternative paths (P4) 
                # or suppose it has alternative paths set b \in R' but 
                # (i) the flow of all ELD paths is 0 
                # (ii) the coefficient of this zero-flow ELI path in any b \in R' is positive 
                # this ELI path and all its alternative paths should be removed to the NUE set
                
                # case 3: replace it with an eligible ELD paths in R'
                # if there is b \in R' s.t. flow > 0, then choose it 
                # otherwise, choose a b \in R' s.t. the coefficient of this zero-flow ELI path in b is negative
                # println("Ra size: ", length(insetRA))
                rA = insetRA[1]
                # for rA in insetRA
                r = insetA[rA]
                indbb = findall(x -> x != 0, coes[rA, :])
                Rp = insetB[indbb]
                # println(coes[rA, indbb])
                cond4 = length(indbb) == 0
                if cond4 
                    insetA = setdiff(insetA, r)
                    insetC = union(insetC, r)
                    pi_path1[r] = 0.001
                    pi_path1[Rp] .= 0.001
                    if verbose > 1; println("case 2") end
                else 
                    cond51 = (sum(flow_path1[Rp] .== 0) == length(indbb))
                    cond52 = (sum(coes[rA, indbb] .< 0) == 0)
                    if cond51 & cond52 
                        insetA = setdiff(insetA, r)
                        insetB = setdiff(insetB, Rp)
                        insetC = union(insetC, r, Rp)
                        pi_path1[r] = 0.001
                        pi_path1[Rp] .= 0.001
                        if verbose > 1; println("case 2") end
                    elseif !cond51 
                        indbb = findall(x -> x != 0, flow_path1[Rp])[1]
                        indb = Rp[indbb]
                        insetA = setdiff(insetA, r)
                        insetA = union(insetA, indb)
                        insetB = setdiff(insetB, indb)
                        insetB = union(insetB, r)
                        if verbose > 1
                        println("case 3a, new flow is ", flow_path1[indb])
                        end
                    elseif !cond52 
                        # case 3
                        indbb = (indbb[findall(x -> x < 0, coes[rA, indbb])])[1]
                        indb = insetB[indbb]
                        insetA = setdiff(insetA, r)
                        insetA = union(insetA, indb)
                        insetB = setdiff(insetB, indb)
                        insetB = union(insetB, r)
                        flow_path1[indb] = 0.001
                        if verbose > 1
                        println("case 3b, new flow is ", flow_path1[indb])
                        end
                    else 
                        @assert true println("edge case!!!!!! cond4 is ", cond4, " cond51 is ", cond51, " cond52 is ", cond52)
                    end
                end
            end

            # check if desent condition is satisfied
            ep1.y = ep0.y - g * alpha 
            Fnew = upper_obj(td, flow1, ep1)
            if verbose > 1
                @printf("new obj is %.2f, old obj is %.2f\n", Fnew, Fold)
            end
        
            if Fnew < Fold
                push!(alphas, alphas[end] + alpha)
                push!(Fs, Fnew)
                Fold = Fnew
                flow_path0 = deepcopy(flow_path1)
                pi_path0 = deepcopy(pi_path1)
                flow0 = deepcopy(flow1)
                ep0 = deepcopy(ep1)
            else
                decent_flag = false
            end
            if length(alphas) >= 20; decent_flag = false; end
            if T + (time() - t0) >= 1000; decent_flag = false; end
        end

        # update the solution
        stop_flag = false
        for alpha in reverse(alphas)
            flow0 = deepcopy(flow)
            m_od0 = deepcopy(m_od)
            ep0 = deepcopy(ep)
            ep0.y -= alpha * g
            ep0.y = max.(ep0.y, 0)
            
            for _ in 1:3
                flow0 = ISMO_f!(m_od0, flow0, td, ep0, graph, link_dic, qls_opt)
            end 
            if verbose > 1
                @printf("alpha = %.2e, new obj is %.2f\n", alpha, upper_obj(td, flow0, ep0))
            end
            if alpha == 0
                stop_flag = true 
                if verbose > 0; println("step size is zero, break"); end 
            end
            if upper_obj(td, flow0, ep0) < upper_obj(td, flow, ep)
                ep = deepcopy(ep0)
                flow = deepcopy(flow0)
                m_od = deepcopy(m_od0)
                if verbose > 0; @printf("success update with alpha: %.2e\n", alpha); end
                break
            end
            
        end
        flow = deepcopy(flow0)
        ep = deepcopy(ep0)
        T += time() - t0
        append!(df0, DataFrame(iter = nIter, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
        push!(ys, deepcopy(ep.y))
        println()
        @printf("==============iteration: %4d upper obj: %.4f time: %.2f =================== \n", nIter, upper_obj(td, flow, ep), T)
        if stop_flag; break; end
        if T >= 1000; break; end
        GC.gc()
    end

    CSV.write(results_root*name*"_EBA_jam_$(jam)_per_$(per).csv", df0)
    CSV.write(results_root*name*"_EBA_y_jam_$(jam)_per_$(per).csv", DataFrame(ys, :auto))

    if verbose > 0
        for y in ep.y
            @printf(" %6.2f ", y)
        end
        println()
    end

end

function main_NRMFD_EBA()
jam = 1.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
NRMFD_EBA(;name = "Anaheim", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2., dstd = 4., verbose=1, due_acc = 1e-5)
end
jam = 3.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
NRMFD_EBA(;name = "Anaheim", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2., dstd = 4., verbose=1, due_acc = 1e-4)
end
jam = 5.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
NRMFD_EBA(;name = "Anaheim", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2., dstd = 4., verbose=1, due_acc = 1e-3)
end

jam = 1.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
NRMFD_EBA(;name = "Barcelona", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2e6, dstd = 4e6, verbose=1, due_acc = 1e-4)
end
jam = 3.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
NRMFD_EBA(;name = "Barcelona", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2e6, dstd = 4e6, verbose=1, due_acc = 1e-2)
end
jam = 5.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
NRMFD_EBA(;name = "Barcelona", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2e6, dstd = 4e6, verbose=1, due_acc = 1e-2)
end

jam = 1.0
for per in [0.1, 0.3, 0.7]
println("jam factor: ", jam, "name", name, "per", per)
NRMFD_EBA(;name = "Chicago-Sketch", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2e2, dstd = 4e2, verbose=1, due_acc = 1e-4)
end
jam = 3.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
println("jam factor: ", jam, "name", name, "per", per)
NRMFD_EBA(;name = "Chicago-Sketch", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2e2, dstd = 4e2, verbose=1, due_acc = 1e-2)
end
jam = 5.0
println("jam factor: ", jam)
for per in [0.1, 0.3, 0.7]
println("jam factor: ", jam, "name", name, "per", per)
NRMFD_EBA(;name = "Chicago-Sketch", jam=jam, maxIter = 10, per = per, seed = 1, dmean = 2e2, dstd = 4e2, verbose=1, due_acc = 1e-2)
end
end

main_NRMFD_EBA()