results_root = "../results/"
due_acc = 1e-10

function SiouxFalls_initial(; jam=1.)
    td, flow, ep = load_SiouxFalls(jam)
    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 2)
    
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end

    all_or_nothing!(m_od, flow, td, graph, link_dic);
    
    ep.y = 0 * ones(10)
    while rgap(flow, td, graph, link_dic, ep) > 1e-10
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
    end
    df0 = DataFrame(method="DC", time = 0., upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    append!(df0, DataFrame(method="EBA", time = 0., upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))

    ep.y = 12.5 * ones(10)
    while rgap(flow, td, graph, link_dic, ep) > 1e-10
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
    end
    append!(df0, DataFrame(method="EDO", time = 0., upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
    append!(df0, DataFrame(method="AL", time = 0., upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))


    ep.y = 2 * ones(10)
    while rgap(flow, td, graph, link_dic, ep) > 1e-10
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
    end
    append!(df0, DataFrame(method="HJ", time = 0., upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))

    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_initial_jam_$(jam).csv", df0)
end

function post_SiouxFall(; jam=1., method = "DC")
    td, flow, ep = load_SiouxFalls(jam)
    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 2)
    
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
    
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    while rgap(flow, td, graph, link_dic, ep) > 1e-10
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
    end
    df0 = DataFrame(time = 0., upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))

    ys = Matrix(CSV.read(results_root*"SiouxFalls_results/SiouxFalls_$(method)_y_jam_$(jam).csv", DataFrame))
    times = CSV.read(results_root*"SiouxFalls_results/SiouxFalls_$(method)_jam_$(jam).csv", DataFrame).time
    n_y = size(ys, 2)

    for i in 1:n_y
        ep.y = ys[:, i]
        while rgap(flow, td, graph, link_dic, ep) > 1e-10
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
        end
        append!(df0, DataFrame(time = times[i], upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
    end

    CSV.write(results_root*"SiouxFalls_results/post_$(method)_jam_$(jam).csv", df0)
end

function comparison_SiouxFall(;jam = 1)
    td, flow, ep = load_SiouxFalls(jam)
    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 2)
    
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
    
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    while rgap(flow, td, graph, link_dic, ep) > 1e-10
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
        # @printf("Upper Obj is %8.6f (rgap = %8.4e)\n", upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
    end
    df0 = DataFrame(method="None", time = 0., upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))

    for method in ["DC", "DCE", "HJ", "EDO", "AL", "EBA"]
        td, flow, ep = load_SiouxFalls(jam)
        graph = create_graph(td.init_node, td.term_node)
        link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))
    
        qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
        
        m_od = Array{Any}(nothing, size(td.travel_demand)[1])
        for i in 1:size(td.travel_demand)[1]
            m_od[i] = Array{Any}(nothing, 0)
        end
        
        ep.y = CSV.read(results_root*"SiouxFalls_results/SiouxFalls_$(method)_y_jam_$(jam).csv", DataFrame)[:, end]
        all_or_nothing!(m_od, flow, td, graph, link_dic);
        while rgap(flow, td, graph, link_dic, ep) > 1e-10
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
            # @printf("Upper Obj is %8.6f (rgap = %8.4e)\n", upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
        end

        append!(df0, DataFrame(method=method, time = CSV.read(results_root*"SiouxFalls_results/SiouxFalls_$(method)_jam_$(jam).csv", DataFrame).time[end], upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
    end

    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_comparison_jam_$(jam).csv", df0)

end


function eval_ep_ISMO(td, m_od, flow, ep, graph, link_dic, qls_opt)
    for _ in 1:10
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
    end
    result = upper_obj(td, flow, ep)
    @printf("Upper Obj is %8.6f", result)
end

# SiouxFalls
function DC_SiouxFalls(;y0 = zeros(10), jam = 1)
    @printf("\n\n here is DC method\n")

    # for precompiling
    td, flow, ep = load_SiouxFalls(jam)
    # y0 = 2. * ones(10)
    ep.y = deepcopy(y0)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 2)
    m_od_f = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od_f[i] = Array{Any}(nothing, 0)
    end

    # Initialization
    T = 0
    timer = time()
    all_or_nothing!(m_od_f, flow, td, graph, link_dic);
    for _ in 1:1
    # while rgap(flow, td, graph, link_dic, ep) > due_acc
        flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
    end
    T += time() - timer 
    ys = [deepcopy(ep.y)]
    
    # for pre compiling
    
    timer = time()
    m_od_F = deepcopy(m_od_f)
    Flow = deepcopy(flow)
    
    T += time() - timer
    df0 = DataFrame(iter = 0, time = T, upper1 = upper_obj(td, Flow, ep), upper2 = upper_obj(td, flow, ep), aec1 = aec(Flow, td, graph, link_dic, ep), aec2 = aec(flow, td, graph, link_dic, ep), xi = 0., err = 0., rho = bi_para.rho)
    @printf("%3d%-24s %.4fs \t upper_obj: %8.6f %8.6f (aec = %9.2e, %9.2e, xi = %9.2e, err = %9.2e)\n", 0, "-Iterations take", T, upper_obj(td, Flow, ep), upper_obj(td, flow, ep), aec(Flow, td, graph, link_dic, ep), aec(flow, td, graph, link_dic, ep), 0., 0.)
    
    for iter in 1:1
        
        timer = time()
        ep0 = deepcopy(ep)
        Flow0 = deepcopy(Flow)

        for _ in 1:1
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
        end
        
        flag = true 
        while flag
            y1 = deepcopy(ep.y)
            Flow = ISMO_F!(m_od_F, Flow, Flow0, td, ep, graph, link_dic, qls_opt, bi_para)
            enhance3!(ep, ep0.y, Flow, flow, td, bi_para)
            flag = (norm(ep.y - y1) > 1e-10) 
        end

        xi = cal_xi_dc(td, Flow, ep, flow, ep0)
        # xi = cal_xi_dc(td, Flow, ep, flow, ep0, m_od_f, graph, link_dic, qls_opt)
        err = err_vy(Flow, ep.y, Flow0, ep0.y)
        T += time() - timer 
        @printf("%3d%-24s %.4fs \t upper_obj: %8.6f %8.6f (aec = %9.2e, %9.2e, xi = %9.2e, err = %9.2e)\n", iter, "-Iterations take", T, upper_obj(td, Flow, ep), upper_obj(td, flow, ep0), aec(Flow, td, graph, link_dic, ep), aec(flow, td, graph, link_dic, ep0), xi, err)
        append!(df0, DataFrame(iter = iter, time = T, upper1 = upper_obj(td, Flow, ep), upper2 = upper_obj(td, flow, ep0), aec1 = aec(Flow, td, graph, link_dic, ep), aec2 = aec(flow, td, graph, link_dic, ep0), xi = xi, err = err, rho = bi_para.rho))
        push!(ys, deepcopy(ep.y))

        # if abs(xi) < 1e-4; break; end
        if (xi) < 1e-3; break; end
        # if df0[end, :upper2] > df0[end-1, :upper2] + 1e-2; break; end

        timer = time()
        if max(bi_para.rho, 1/xi) < 1/err
            bi_para.rho *= bi_para.tau  
            @printf("rho update to %4f\n", bi_para.rho)
        end
        T += time() - timer 
    end

    # for real_calculation
    td, flow, ep = load_SiouxFalls(jam)
    ep.y = deepcopy(y0)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 2)
    m_od_f = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od_f[i] = Array{Any}(nothing, 0)
    end

    # Initialization
    T = 0
    timer = time()
    all_or_nothing!(m_od_f, flow, td, graph, link_dic);
    # for _ in 1:10
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
    end
    T += time() - timer 
    ys = [deepcopy(ep.y)]
    
    # for pre compiling
    
    timer = time()
    m_od_F = deepcopy(m_od_f)
    Flow = deepcopy(flow)
    
    T += time() - timer
    df0 = DataFrame(iter = 0, time = T, upper1 = upper_obj(td, Flow, ep), upper2 = upper_obj(td, flow, ep), aec1 = aec(Flow, td, graph, link_dic, ep), aec2 = aec(flow, td, graph, link_dic, ep), xi = 0., xi_exact = 0., err = 0., rho = bi_para.rho)

    @printf("%3d%-24s %.4fs \t upper_obj: %8.6f %8.6f (aec = %9.2e, %9.2e, xi = %9.2e, err = %9.2e)\n", 0, "-Iterations take", T, upper_obj(td, Flow, ep), upper_obj(td, flow, ep), aec(Flow, td, graph, link_dic, ep), aec(flow, td, graph, link_dic, ep), 0., 0.)
    
    for iter in 1:200
        
        timer = time()
        ep0 = deepcopy(ep)
        Flow0 = deepcopy(Flow)

        while rgap(flow, td, graph, link_dic, ep) > due_acc
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt, false)
        end
        
        flag = true 
        while flag
            y1 = deepcopy(ep.y)
            Flow = ISMO_F!(m_od_F, Flow, Flow0, td, ep, graph, link_dic, qls_opt, bi_para)
            enhance3!(ep, ep0.y, Flow, flow, td, bi_para)
            flag = (norm(ep.y - y1) > 1e-5) 
        end

        xi = cal_xi_dc(td, Flow, ep, flow, ep0)
        err = err_vy(Flow, ep.y, Flow0, ep0.y)
        T += time() - timer 

        ep1 = deepcopy(ep)
        Flow1 = deepcopy(Flow)
        m_od_F1 = deepcopy(m_od_F)
        while rgap(Flow1, td, graph, link_dic, ep1) > due_acc
            Flow1 = ISMO_f!(m_od_F1, Flow1, td, ep1, graph, link_dic, qls_opt)
            Flow1 = ISMO_f!(m_od_F1, Flow1, td, ep1, graph, link_dic, qls_opt, false)
        end        
        fk1 = lower_obj(td, Flow, ep1)
        vk1 = lower_obj(td, Flow1, ep1)
        xi_exact = upper_obj(td, Flow, ep1) + bi_para.rho * (fk1 - vk1)

        @printf("%3d%-24s %.4fs \t upper_obj: %8.6f %8.6f (aec = %9.2e, %9.2e, xi = %9.2e, err = %9.2e)\n", iter, "-Iterations take", T, upper_obj(td, Flow, ep), upper_obj(td, flow, ep0), aec(Flow, td, graph, link_dic, ep), aec(flow, td, graph, link_dic, ep0), xi, err)
        append!(df0, DataFrame(iter = iter, time = T, upper1 = upper_obj(td, Flow, ep), upper2 = upper_obj(td, flow, ep0), aec1 = aec(Flow, td, graph, link_dic, ep), aec2 = aec(flow, td, graph, link_dic, ep0), xi = xi, xi_exact = xi_exact, err = err, rho = bi_para.rho))
        push!(ys, deepcopy(ep.y))

        # if abs(xi) < 1e-4; break; end
        if (xi) < 1e-3; break; end

        timer = time()
        if max(bi_para.rho, 1/xi) < 1/err
            bi_para.rho *= bi_para.tau
            @printf("rho update to %4f\n", bi_para.rho)
        end
        T += time() - timer 
    end
    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_DC_jam_$(jam).csv", df0)
    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_DC_y_jam_$(jam).csv", DataFrame(ys, :auto))

    @printf("%24s %.4fs \t upper_obj: %8.6f (aec = %6.2e)\n", "Iteration takes", T, upper_obj(td, Flow, ep), aec(Flow, td, graph, link_dic, ep))
    for i in 1:lastindex(ep.y); @printf("%6.2f, ", ep.y[i]); end
    @printf("\n")
    for i in 1:lastindex(ep.y); @printf("%6.2f& ", ep.y[i]); end
    @printf("\\\\ \n")

    Eval_EP_SiouxFall(ep.y, jam = jam)
end


function EDO_SiouxFalls(;jam=1)
    @printf("\n\n here is EDO method\n")
    # for precompiling
    td, flow, ep = load_SiouxFalls(jam)

    EPS_EDO = 1e-5
    @views inds = ep.inds
    ep.y = 12.5 * ones(lastindex(inds))
    yL = zeros(lastindex(inds))
    yU = 25 * ones(lastindex(inds))

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
 
    # Initialization
    T = 0
    timer = time()
    
    # EDO 
    @views inds = ep.inds
    @views d = ep.d
    yy = (yL + yU) / 2
    ep.y = yy
    iter = 0
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    for _ in 1:1
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
    end
    
    T += time() - timer 
    @printf("%5d-th iteration: %6.2fs upper_obj: %8.6f (aec = %5.2e)\n", 0, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
    df0 = DataFrame(iter = 0, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    ys = [deepcopy(yy)]
    
    flag = true
    for i in 1:1
        timer = time()
        flag = false
        iter += 1
        for i in 1:lastindex(inds)
            if ~flag 
                if yU[i] - yL[i] < EPS_EDO; 
                    continue
                else
                    flag = true
                end
            end

            yy[i] = (yL[i]+yU[i])/2
            tmp = EDO_dZa(td, inds[i], d[i], flow[inds[i]], yy[i])

            if tmp < 0
                yL[i] = yy[i] 
            elseif tmp > 0 
                yU[i] = yy[i] 
            else
                yL[i] = yy[i]
                yU[i] = yy[i]
            end
        end
        ep.y = yy
        for _ in 1:1
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        end
        T += time() - timer
        @printf("%5d-th iteration: %6.2fs upper_obj: %8.6f (aec = %5.2e)\n", iter, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
        append!(df0, DataFrame(iter = iter, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
        push!(ys, deepcopy(yy))
    end

    # for real calculation
    td, flow, ep = load_SiouxFalls(jam)

    EPS_EDO = 1e-5
    @views inds = ep.inds
    ep.y = 12.5 * ones(lastindex(inds))
    yL = zeros(lastindex(inds))
    yU = 25 * ones(lastindex(inds))

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
 
    # Initialization
    T = 0
    timer = time()
    
    # EDO 
    @views inds = ep.inds
    @views d = ep.d
    yy = (yL + yU) / 2
    ep.y = yy
    iter = 0
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
    end
    
    T += time() - timer 
    @printf("%5d-th iteration: %6.2fs upper_obj: %8.6f (aec = %5.2e)\n", 0, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
    df0 = DataFrame(iter = 0, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    ys = [deepcopy(yy)]
    
    flag = true
    while flag 
        timer = time()
        flag = false
        iter += 1
        for i in 1:lastindex(inds)
            if ~flag 
                if yU[i] - yL[i] < EPS_EDO; 
                    continue
                else
                    flag = true
                end
            end

            yy[i] = (yL[i]+yU[i])/2
            tmp = EDO_dZa(td, inds[i], d[i], flow[inds[i]], yy[i])

            if tmp < 0
                yL[i] = yy[i] 
            elseif tmp > 0 
                yU[i] = yy[i] 
            else
                yL[i] = yy[i]
                yU[i] = yy[i]
            end
        end
        ep.y = yy
        # for _ in 1:3
        while rgap(flow, td, graph, link_dic, ep) > due_acc
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
        end
        T += time() - timer
        @printf("%5d-th iteration: %6.2fs upper_obj: %8.6f (aec = %5.2e)\n", iter, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
        append!(df0, DataFrame(iter = iter, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
        push!(ys, deepcopy(yy))
    end

    ep.y = (yL + yU) / 2
    @printf("%24s %.4fs \t upper_obj: %8.6f (aec = %6.2e)\n", "Iteration takes", T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
    append!(df0, DataFrame(iter = iter+1, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_EDO_jam_$(jam).csv", df0)
    
    push!(ys, deepcopy(ep.y))
    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_EDO_y_jam_$(jam).csv", DataFrame(ys, :auto))

    @printf("%24s %.4fs \t upper_obj: %8.6f (aec = %6.2e)\n", "Iteration takes", T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
    for i in 1:lastindex(ep.y); @printf("%6.2f, ", ep.y[i]); end
    @printf("\n")
    for i in 1:lastindex(ep.y); @printf("%6.2f& ", ep.y[i]); end
    @printf("\\\\ \n")

    Eval_EP_SiouxFall(ep.y, jam = jam)
end

function HJ_SiouxFalls(;y0 = 2*ones(10), jam = 1)
    @printf("\n\n here is HJ method\n")
    # for pre compiling
    td, flow, ep = load_SiouxFalls(jam)
    para = MutableNamedTuple(eps = 1e-1, alpha = 0.9, delta = .2)

    ep.y = deepcopy(y0)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
 
    # Initialization
    T = 0 
    timer = time()
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    for _ in 1:1
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
    end
    T += time() - timer
    df0 = DataFrame(iter = 0, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    ys = [deepcopy(ep.y)]

    # HJ method
    timer = time()
    iter = 0
    z = deepcopy(ep.y)
    @views inds = ep.inds
    @views d = ep.d
    T += time() - timer
    @printf("%5d-th iteration: %6.4fs upper_obj: %8.6f (aec = %5.2e)\n", iter, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))

    # while true
    for i in 1:1
        timer = time()
        iter += 1 
        flow_z = deepcopy(flow)
        m_od_z = deepcopy(m_od)
        ep_z = deepcopy(ep)
        for i in 1:lastindex(inds)
            if HJ_Fz(td, flow_z, inds[i], d[i], z[i] + para.delta ) < HJ_Fz(td, flow_z, inds[i], d[i], z[i])
                z[i] = z[i] + para.delta
                ep_z.y = z
            # while rgap(flow_z, td, graph, link_dic, ep) > due_acc
            # for i in 1:1
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt)
            # end
            elseif HJ_Fz(td, flow_z, inds[i], d[i], z[i] - para.delta ) < HJ_Fz(td, flow_z, inds[i], d[i], z[i])
                z[i] = z[i] - para.delta
                ep_z.y = z
            # while rgap(flow_z, td, graph, link_dic, ep) > due_acc
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt)
            # end
            end 
        end

        if HJ_Fy(td, flow_z, ep_z) < HJ_Fy(td, flow, ep)
            y0 = deepcopy(ep.y)
            ep.y = deepcopy(z)
            for _ in 1:1
               flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
            end
            z = deepcopy(ep.y) + para.alpha * (ep.y - y0)
        elseif para.delta <= para.eps
            T += time() - timer 
            break 
        else 
            para.delta *= 0.5
            z = deepcopy(ep.y)
        end

        T += time() - timer 
        @printf("%5d-th iteration: %6.4fs upper_obj: %8.6f (aec = %5.2e)\n", iter, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
        append!(df0, DataFrame(iter = iter, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
        push!(ys, deepcopy(ep.y))

    end

    # for real calculation 
    td, flow, ep = load_SiouxFalls(jam)
    para = MutableNamedTuple(eps = 1e-1, alpha = 0.9, delta = .2)

    ep.y = deepcopy(y0)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
 
    # Initialization
    T = 0 
    timer = time()
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    # for _ in 1:10
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
    end
    T += time() - timer
    df0 = DataFrame(iter = 0, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    ys = [deepcopy(ep.y)]

    # HJ method
    timer = time()
    iter = 0
    z = deepcopy(ep.y)
    @views inds = ep.inds
    @views d = ep.d
    T += time() - timer
    @printf("%5d-th iteration: %6.4fs upper_obj: %8.6f (aec = %5.2e)\n", iter, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))

    while true
        timer = time()
        iter += 1 
        flow_z = deepcopy(flow)
        m_od_z = deepcopy(m_od)
        ep_z = deepcopy(ep)
        for i in 1:lastindex(inds)
            if HJ_Fz(td, flow_z, inds[i], d[i], z[i] + para.delta ) < HJ_Fz(td, flow_z, inds[i], d[i], z[i])
                z[i] = z[i] + para.delta
                ep_z.y = z
            while rgap(flow_z, td, graph, link_dic, ep) > due_acc
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt)
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt, false)
            end
            elseif HJ_Fz(td, flow_z, inds[i], d[i], z[i] - para.delta ) < HJ_Fz(td, flow_z, inds[i], d[i], z[i])
                z[i] = z[i] - para.delta
                ep_z.y = z
            while rgap(flow_z, td, graph, link_dic, ep) > due_acc
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt)
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt, false)
            end
            end 
        end

        if HJ_Fy(td, flow_z, ep_z) < HJ_Fy(td, flow, ep)
            y0 = deepcopy(ep.y)
            ep.y = deepcopy(z)
            # for _ in 1:1
            while rgap(flow, td, graph, link_dic, ep) > due_acc
               flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
               flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
            end
            z = deepcopy(ep.y) + para.alpha * (ep.y - y0)
        elseif para.delta <= para.eps
            T += time() - timer 
            break 
        else 
            para.delta *= 0.5
            z = deepcopy(ep.y)
        end

        T += time() - timer 
        @printf("%5d-th iteration: %6.4fs upper_obj: %8.6f (aec = %5.2e)\n", iter, T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
        append!(df0, DataFrame(iter = iter, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
        push!(ys, deepcopy(ep.y))

    end

    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_HJ_jam_$(jam).csv", df0)
    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_HJ_y_jam_$(jam).csv", DataFrame(ys, :auto))

    @printf("%24s %.4fs \t upper_obj: %8.6f (aec = %6.2e)\n", "Iteration takes", T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
    for i in 1:lastindex(ep.y); @printf("%6.2f, ", ep.y[i]); end
    @printf("\n")
    for i in 1:lastindex(ep.y); @printf("%6.2f& ", ep.y[i]); end
    @printf("\\\\ \n")

    Eval_EP_SiouxFall(ep.y, jam = jam)
end

function AL_SiouxFalls(;jam = 1)
    @printf("\n\n here is AL method\n")
    # for pre compiling
    td, flow, ep = load_SiouxFalls(jam)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    # m_od = Array{Any}(nothing, size(td.travel_demand)[1], size(td.travel_demand)[2])
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
 
    # Initialization
    al_para = MutableNamedTuple(mu = 0., rho = 1, beta = 100, EPS_AL = 1e-2, EPS_IN = 5e-2)
    T = 0
    # AL 
    @views inds = ep.inds
    ep.y = 12.5 * ones(lastindex(inds))
    yL = zeros(lastindex(inds))
    yU = 25 * ones(lastindex(inds))

    timer = time()
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    for i = 1:1
    # while rgap(flow, td, graph, link_dic, ep) > due_acc
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
    end
    T += time() - timer

    df0 = DataFrame(iter = 0, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    ys = [deepcopy(ep.y)]

    timer = time()
    flow2 = deepcopy(flow)
    h0 = lower_obj(td, flow, ep) - lower_obj(td, flow2, ep)
    rho = al_para.rho
    mu = al_para.mu
    T += time() - timer 

    k = 0
    h2 = 0
    for i = 1:1
        timer = time()
        # 1 inner loop 
        mu = al_para.mu
        rho = al_para.rho
        EPS_IN = al_para.EPS_IN

        k = 0
        for i = 1:1
            k = k+1 
            # 1 DUE 
            for i in 1:1
                flow2 = ISMO_f!(m_od, flow2, td, ep, graph, link_dic, qls_opt)
                flow2 = ISMO_f!(m_od, flow2, td, ep, graph, link_dic, qls_opt, false)
            end

        
            # 2 calculation marginal    
            # 3 calculation general cost 
            tk, gk, h2 = AL_tg(flow, flow2, ep, td, mu, rho)

            @assert minimum(tk) >= 0 "the minimum is $(minimum(tk))"
            
            # 4 all or nothing for v / y 
            flow3 = all_or_nothing(tk, td, graph, link_dic)
            y3 = zeros(lastindex(ep.inds))
            for i in 1:lastindex(ep.inds)
                if gk[i] > 0; y3[i] = yL[i]
                else; y3[i] = yU[i]
                end
            end

            # 5 stopping test
            # if k > 200; break; end
            d1 = flow3 - flow 
            d2 = y3 - ep.y 
            Zk = abs(tk'*d1 + gk'*d2) / (upper_obj(td, flow, ep) + mu*h2 + .5*rho*h2^2)
            # println(Zk)
            if Zk < EPS_IN; break; end 

            # 6 line search
            alpha = 1 / (k + 10)

            # 7 update 
            flow = flow + alpha * d1 
            ep.y = ep.y + alpha * d2
            

        end 
        # return flow, ep, h2

        # 2 stopping test
        if h2 < al_para.EPS_AL
            T += time() - timer
            @printf("%5d-th iteration takes %6.2fs upper_obj: %8.6f (aec = %5.2e, h = %5.2e)\n", i , T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep), h2)
            append!(df0, DataFrame(iter = i, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
            push!(ys, deepcopy(ep.y))
            break
        end

        # 3 updating Lagrangian multiplie
        # 4 update penalty parameter
        if h2 < 0.25 * h0
            al_para.mu = al_para.mu + al_para.rho * h2
        else
            al_para.rho = al_para.beta * al_para.rho 
        end
        h0 = h2

        T += time() - timer
        println("mu = ", al_para.mu, "  rho = ", al_para.rho)
        @printf("%5d-th iteration takes %6.2fs upper_obj: %8.6f (aec = %5.2e, h = %5.2e)\n", i , T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep), h2)
        append!(df0, DataFrame(iter = i, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
        push!(ys, deepcopy(ep.y))
    end

    # for real calculation 
    td, flow, ep = load_SiouxFalls(jam)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    # m_od = Array{Any}(nothing, size(td.travel_demand)[1], size(td.travel_demand)[2])
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
 
    # Initialization
    al_para = MutableNamedTuple(mu = 0., rho = 1, beta = 100, EPS_AL = 1e-2, EPS_IN = 5e-2)
    T = 0
    # AL 
    @views inds = ep.inds
    ep.y = 12.5 * ones(lastindex(inds))
    yL = zeros(lastindex(inds))
    yU = 25 * ones(lastindex(inds))

    timer = time()
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
    end
    T += time() - timer

    df0 = DataFrame(iter = 0, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep))
    ys = [deepcopy(ep.y)]

    timer = time()
    flow2 = deepcopy(flow)
    h0 = lower_obj(td, flow, ep) - lower_obj(td, flow2, ep)
    rho = al_para.rho
    mu = al_para.mu
    T += time() - timer 

    k = 0
    h2 = 0
    for i = 1:10
        timer = time()
        # 1 inner loop 
        mu = al_para.mu
        rho = al_para.rho
        EPS_IN = al_para.EPS_IN

        k = 0
        while true
            k = k+1 
            # 1 DUE 
            while rgap(flow2, td, graph, link_dic, ep) > 1e-10
                flow2 = ISMO_f!(m_od, flow2, td, ep, graph, link_dic, qls_opt)
                flow2 = ISMO_f!(m_od, flow2, td, ep, graph, link_dic, qls_opt, false)
            end

        
            # 2 calculation marginal    
            # 3 calculation general cost 
            tk, gk, h2 = AL_tg(flow, flow2, ep, td, mu, rho)

            @assert minimum(tk) >= 0 "the minimum is $(minimum(tk))"
            
            # 4 all or nothing for v / y 
            flow3 = all_or_nothing(tk, td, graph, link_dic)
            y3 = zeros(lastindex(ep.inds))
            for i in 1:lastindex(ep.inds)
                if gk[i] > 0; y3[i] = yL[i]
                else; y3[i] = yU[i]
                end
            end

            # 5 stopping test
            # if k > 200; break; end
            d1 = flow3 - flow 
            d2 = y3 - ep.y 
            Zk = abs(tk'*d1 + gk'*d2) / (upper_obj(td, flow, ep) + mu*h2 + .5*rho*h2^2)
            # println(Zk)
            if Zk < EPS_IN; break; end 

            # 6 line search
            alpha = 1 / (k + 10)

            # 7 update 
            flow = flow + alpha * d1 
            ep.y = ep.y + alpha * d2
            
        end 
        # return flow, ep, h2

        # 2 stopping test
        if h2 < al_para.EPS_AL
            T += time() - timer
            @printf("%5d-th iteration takes %6.2fs upper_obj: %8.6f (aec = %5.2e, h = %5.2e)\n", i , T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep), h2)
            append!(df0, DataFrame(iter = i, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
            push!(ys, deepcopy(ep.y))
            break
        end

        # 3 updating Lagrangian multiplie
        # 4 update penalty parameter
        if h2 < 0.25 * h0
            al_para.mu = al_para.mu + al_para.rho * h2
        else
            al_para.rho = al_para.beta * al_para.rho 
        end
        h0 = h2

        T += time() - timer
        println("mu = ", al_para.mu, "  rho = ", al_para.rho)
        @printf("%5d-th iteration takes %6.2fs upper_obj: %8.6f (aec = %5.2e, h = %5.2e)\n", i , T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep), h2)
        append!(df0, DataFrame(iter = i, time = T, upper = upper_obj(td, flow, ep), aec = aec(flow, td, graph, link_dic, ep)))
        push!(ys, deepcopy(ep.y))
    end
    @printf("%24s %.2fs upper_obj: %8.6f (aec = %5.2e)\n", "Iteration takes", T, upper_obj(td, flow, ep), aec(flow, td, graph, link_dic, ep))
    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_AL_jam_$(jam).csv", df0)
    CSV.write(results_root*"SiouxFalls_results/SiouxFalls_AL_y_jam_$(jam).csv", DataFrame(ys, :auto))

    for i in 1:lastindex(ep.y); @printf("%6.2f, ", ep.y[i]); end
    @printf("\n")
    for i in 1:lastindex(ep.y); @printf("%6.2f& ", ep.y[i]); end
    @printf("\\\\ \n")

    Eval_EP_SiouxFall(ep.y, jam = jam)
end


## DC with data saved
function DC_SiouxFalls_data(;y0 = zeros(10))
    name = "SiouxFalls"
    open("results/DC_" * name * ".txt", "w") do io; end
    open("results/DC_" * name * "_y.txt", "w") do io; end
    open("results/DC_" * name * "_flowU.txt", "w") do io; end
    open("results/DC_" * name * "_flowL.txt", "w") do io; end
    # this part is for precompiling 
    td, flow, ep = load_SiouxFalls()
    # y0 = 2. * ones(10)
    ep.y = deepcopy(y0)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 100)
    # m_od_f = Array{Any}(nothing, size(td.travel_demand)[1], size(td.travel_demand)[2])
    m_od_f = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od_f[i] = Array{Any}(nothing, 0)
    end

    # Initialization
    all_or_nothing!(m_od_f, flow, td, graph, link_dic);
    for i in 1:10
        flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
    end
    
    # for pre compiling
    m_od_F = deepcopy(m_od_f)
    Flow = deepcopy(flow)
    
    timer = time()
    time_print = 0

    for iter in 1:1
        
        ep0 = deepcopy(ep)
        Flow0 = deepcopy(Flow)

        # for _ in 1:2
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt, false)
        # end
        
        for _ in 1:1
            enhance3!(ep, ep0.y, Flow, flow, td, bi_para)
            
            Flow = ISMO_F!(m_od_F, Flow, Flow0, td, ep, graph, link_dic, qls_opt, bi_para)
            enhance3!(ep, ep0.y, Flow, flow, td, bi_para)
        end

        xi = cal_xi_dc(td, Flow, ep, flow, ep0)
        err = err_vy(Flow, ep.y, Flow0, ep0.y)
        timep = @elapsed begin
            @printf("%3d%-24s %.4fs \t upper_obj: %8.6f %8.6f (aec = %9.2e, %9.2e, xi = %9.2e, err = %9.2e)\n", iter, "-Iterations take", time() - timer - time_print, upper_obj(td, Flow, ep), upper_obj(td, flow, ep0), aec(Flow, td, graph, link_dic, ep), aec(flow, td, graph, link_dic, ep0), xi, err)
            open("results/DC_" * name * ".txt", "a") do io
                @printf(io, "%8.2f %9.2e %6d\n", time() - timer - time_print, xi, bi_para.rho)
            end
            open("results/DC_" * name * "_y.txt", "a") do io
                for i in 1:lastindex(ep.y); @printf(io, "%6.4f ", ep.y[i]); end
                @printf(io, "\n")
            end
            open("results/DC_" * name * "_flowU.txt", "a") do io
                for i in 1:lastindex(Flow); @printf(io, "%6.4f ", Flow[i]); end
                @printf(io, "\n")
            end
            open("results/DC_" * name * "_flowL.txt", "a") do io
                for i in 1:lastindex(flow); @printf(io, "%6.4f ", flow[i]); end
                @printf(io, "\n")
            end
        end 
        time_print += timep

        if xi < 1e-3; break; end
        if max(bi_para.rho, 1/xi) < 1/err
            bi_para.rho *= bi_para.tau  
        end
    end

    # this part is for comparison
    open("results/DC_" * name * ".txt", "w") do io; end
    open("results/DC_" * name * "_y.txt", "w") do io; end
    open("results/DC_" * name * "_flowU.txt", "w") do io; end
    open("results/DC_" * name * "_flowL.txt", "w") do io; end
    td, flow, ep = load_SiouxFalls()
    ep.y = deepcopy(y0)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 100)
    m_od_f = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od_f[i] = Array{Any}(nothing, 0)
    end

    # Initialization
    timer = time()
    all_or_nothing!(m_od_f, flow, td, graph, link_dic);
    for i in 1:10
        flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
    end
    
    # for pre compiling
    m_od_F = deepcopy(m_od_f)
    Flow = deepcopy(flow)
    
    time_print = 0

    for iter in 1:100 
        
        ep0 = deepcopy(ep)
        Flow0 = deepcopy(Flow)

        # for _ in 1:2
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt, false)
        # end
        
        for _ in 1:1
            enhance3!(ep, ep0.y, Flow, flow, td, bi_para)
            
            Flow = ISMO_F!(m_od_F, Flow, Flow0, td, ep, graph, link_dic, qls_opt, bi_para)
        end

        xi = cal_xi_dc(td, Flow, ep, flow, ep0)
        err = err_vy(Flow, ep.y, Flow0, ep0.y)
        timep = @elapsed begin
            open("results/DC_" * name * ".txt", "a") do io
                @printf(io, "%8.2f %9.2e %6d\n", time() - timer - time_print, xi, bi_para.rho)
            end
            open("results/DC_" * name * "_y.txt", "a") do io
                for i in 1:lastindex(ep.y); @printf(io, "%6.4f ", ep.y[i]); end
                @printf(io, "\n")
            end
            open("results/DC_" * name * "_flowU.txt", "a") do io
                for i in 1:lastindex(Flow); @printf(io, "%6.4f ", Flow[i]); end
                @printf(io, "\n")
            end
            open("results/DC_" * name * "_flowL.txt", "a") do io
                for i in 1:lastindex(flow); @printf(io, "%6.4f ", flow[i]); end
                @printf(io, "\n")
            end
            @printf("%8.2f %8.2f %9.2e %6d\n", time() - timer, time_print, xi, bi_para.rho)
        end 
        time_print += timep
        if max(bi_para.rho, 1/xi) < 1/err
            bi_para.rho += bi_para.tau
        end
    end

end

function SiouxFalls_post_y(;y0 = zeros(10))
    ydata = readdlm("results/DC_SiouxFalls_y.txt")
    Flows = readdlm("results/DC_SiouxFalls_flowU.txt")
    flows = readdlm("results/DC_SiouxFalls_flowL.txt")

    open("results/DC_SiouxFalls_F.txt", "w") do io; end
    open("results/DC_SiouxFalls_FU.txt", "w") do io; end
    open("results/DC_SiouxFalls_FL.txt", "w") do io; end

    eps = 1e-9
    td, flow, ep = load_SiouxFalls()
    ep.y = y0
    ep0 = deepcopy(ep)
    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))
    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
    all_or_nothing!(m_od, flow, td, graph, link_dic);

    for i in 1:100

        y = ydata[i, :]
        ep.y = y 
        flow1 = Flows[i, :]
        flow2 = flows[i, :]

        # Initialization
        while aec(flow, td, graph, link_dic, ep) > eps
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        end

        open("results/DC_SiouxFalls_F.txt", "a") do io
            @printf(io, "%12.6f\n", upper_obj(td, flow, ep))
        end
        open("results/DC_SiouxFalls_FU.txt", "a") do io
            @printf(io, "%12.6f\n", upper_obj(td, flow1, ep))
        end
        open("results/DC_SiouxFalls_FL.txt", "a") do io
            @printf(io, "%12.6f\n", upper_obj(td, flow2, ep0))
        end
        @printf("%3d-iteration: %12.6f %12.6f %12.6f\n", i, upper_obj(td, flow, ep), upper_obj(td, flow1, ep), upper_obj(td, flow2, ep))
        ep0 = deepcopy(ep)
    end

end