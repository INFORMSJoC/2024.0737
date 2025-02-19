# ==============================================================================
# This script implements a set of optimization algorithms that may be used in CNDP.
#
# Functions in this file:
# - ta_fw: Frank-Wolfe optimization algorithm for User-Equilibrium Problem.
# - ISMO_S!: ISMO-based optimizer for Phi
# - ISMO_f!: ISMO-based optimizer for User-Equilibrium Problem
# - ISMO_F!: ISMO-based optimizer for System-Equilibrium Problem
# - enhance3! and enhance_S!: Enhancement methods for updating the traffic flow plan.
# - AL_tg and AL_INNER: Augmented Lagrangian method to iteratively adjust flow distributions and solve for optimal solutions.
#
# Author: Haian Yin
# Date: 2023-07-21  # Date when the code was written or last modified
# ==============================================================================

# optimizer
function ta_fw(td, flow, ep, graph, link_dic)
    # find user optimal by frank-wolfe
    travel_time = bpr(td, flow, ep)
    flow2 = all_or_nothing(travel_time, td, graph, link_dic)
    d = flow2 - flow 
    optk = optimize(tau -> lower_obj(td, flow+tau*d, ep), 0.0, 1.0, GoldenSection())
    tauk = optk.minimizer
    # println(tauk)
    flow = flow + tauk * d 

    @assert minimum(flow) >= 0

    return flow 
end

function ISMO_S!(m_od, flow::Vector{Float64}, td::TA_Data, ep::Enhancement_plan, graph, link_dic, qls_opt, require_generation = true)
    # find system optimal by ISMO 
    if require_generation
        local state::Graphs.DijkstraState{Float64, Int}
    end
    local flow_cost::Vector{Float64}
    n_path_add = 0
    n_path_del = 0
    
    for r in 1:size(td.travel_demand)[1]
        # for each origin node r, find shortest paths to all destination nodes
        flow_cost = dsys(td, flow, ep)

        if require_generation
            # best solution here is bellman-ford method 
            state = TA_dijkstra_shortest_paths(graph, flow_cost, r, td.init_node, td.term_node, td.first_thru_node)
        end 

        # for s in 1:size(td.travel_demand)[2]
        for ss in 1:length(m_od[r])
            s = m_od[r][ss].D
            # pass if no demand
            if td.travel_demand[r, s] == 0
                continue
            end
            
            # add path if needed 
            # @views tmp = m_od[r, s]
            @views tmp = m_od[r][ss]

            if require_generation
                new_path = get_path(state, r, s, link_dic)
                if !(new_path in tmp.paths)
                    tmp.npath += 1
                    push!(tmp.paths, get_path(state, r, s, link_dic))
                    push!(tmp.flows, 0)
                    n_path_add += 1
                end
            end

            # pass if there is only 1 path for the OD pair
            if tmp.npath == 1; continue; end
            # calculate the direction
            ind_i, ind_j = 1, 1
            cost_p = dsys_p(td, tmp.paths[1], flow, ep)
            cost_i = cost_p 
            cost_j = cost_p
            for k in 2:tmp.npath
                cost_p = dsys_p(td, tmp.paths[k], flow, ep)
                if cost_p < cost_i 
                    ind_i = k 
                    cost_i = cost_p
                end
                if cost_p > cost_j 
                    ind_j = k 
                    cost_j = cost_p 
                end
            end
            if ind_i == ind_j; continue; end 

            d_flow = zeros(length(td.init_node))
            d_flow[tmp.paths[ind_i]] .+= 1
            d_flow[tmp.paths[ind_j]] .-= 1
            @views beta_r = tmp.flows[ind_j]

            path_ij = setdiff(union(tmp.paths[ind_i], tmp.paths[ind_j]), intersect(tmp.paths[ind_i], tmp.paths[ind_j]))
            
            # QLS
            fxk = sys_p(flow, td, path_ij, ep)
            
            # alpha = min(beta_r, qls_opt.lambda)
            # alpha = min(beta_r, (cost_j - cost_i)/lower_obj_p_dHd(flow, td, path_ij, ep))
            alpha = beta_r 
            if alpha > eps()
                fxkp = sys_p(flow + alpha * d_flow, td, path_ij, ep)
                while (abs(fxk - fxkp) > eps()) & (fxkp > fxk - qls_opt.gamma * alpha * alpha)
                    # if fxkp < fxk - qls_opt.gamma * alpha * d_flow2
                    alpha = alpha * qls_opt.delta
                    fxkp = sys_p(flow + alpha * d_flow, td, path_ij, ep)
                end
            end

            # flow update and check if there are 0 flow path (which should be deleted)
            tmp.flows[ind_i] += alpha
            tmp.flows[ind_j] -= alpha
            # tmp.flows[ind_j] = max(tmp.flows[ind_j], 0)
            @assert minimum(tmp.flows[ind_j]) >= 0

            if tmp.flows[ind_i] == 0
                tmp.npath -= 1
                deleteat!(tmp.paths, ind_i)
                deleteat!(tmp.flows, ind_i)
                n_path_add -= 1
            elseif tmp.flows[ind_j] == 0
                tmp.npath -= 1
                deleteat!(tmp.paths, ind_j)
                deleteat!(tmp.flows, ind_j)
                n_path_del += 1
            end
            m_od[r][ss] = tmp
            # m_od[r,s] = tmp

            flow += alpha * d_flow
            flow[path_ij[findall(x -> (x<0), flow[path_ij])]] .= 0
        end

    end

    return flow 
end

function ISMO_f!(m_od, flow::Vector{Float64}, td::TA_Data, ep::Enhancement_plan, graph, link_dic, qls_opt, require_generation = true)
    if require_generation
        local state::Graphs.DijkstraState{Float64, Int}
        # local state::Graphs.BellmanFordState{Float64, Int}
    end
    local flow_cost::Vector{Float64}
    n_path_add = 0
    n_path_del = 0
    
    for r in 1:size(td.travel_demand)[1]
        # for each origin node r, find shortest paths to all destination nodes
        flow_cost = bpr(td, flow, ep)
        
        if require_generation
            state = TA_dijkstra_shortest_paths(graph, flow_cost, r, td.init_node, td.term_node, td.first_thru_node)
            # state = TA_bf_shortest_paths(graph, flow_cost, r, td.init_node, td.term_node, td.first_thru_node)
        end 

        for ss in 1:length(m_od[r])
            s = m_od[r][ss].D 
            # pass if no demand
            if td.travel_demand[r, s] == 0
                continue
            end
            
            # add path if needed 
            @views tmp = m_od[r][ss]

            if require_generation
                new_path = get_path(state, r, s, link_dic)
                if !(new_path in tmp.paths)
                    tmp.npath += 1
                    push!(tmp.paths, get_path(state, r, s, link_dic))
                    push!(tmp.flows, 0)
                    n_path_add += 1
                end
            end

            # pass if there is only 1 path for the OD pair
            if tmp.npath == 1; continue; end
            # calculate the direction
            ind_i, ind_j = 1, 1
            cost_p = bpr_p(td, tmp.paths[1], flow, ep)
            cost_i = cost_p 
            cost_j = cost_p
            for k in 2:tmp.npath
                cost_p = bpr_p(td, tmp.paths[k], flow, ep)
                if cost_p < cost_i 
                    ind_i = k 
                    cost_i = cost_p
                end
                if cost_p > cost_j 
                    ind_j = k 
                    cost_j = cost_p 
                end
            end
            if ind_i == ind_j; continue; end 

            d_flow = zeros(length(td.init_node))
            d_flow[tmp.paths[ind_i]] .+= 1
            d_flow[tmp.paths[ind_j]] .-= 1
            @views beta_r = tmp.flows[ind_j]

            path_ij = setdiff(union(tmp.paths[ind_i], tmp.paths[ind_j]), intersect(tmp.paths[ind_i], tmp.paths[ind_j]))
            
            # QLS
            fxk = lower_obj_p(flow, td, path_ij, ep)
            
            # alpha = min(beta_r, qls_opt.lambda)
            alpha = min(beta_r, (cost_j - cost_i)/lower_obj_p_dHd(flow, td, path_ij, ep))
            # alpha = beta_r
            if alpha > eps()
                fxkp = lower_obj_p(flow + alpha * d_flow, td, path_ij, ep)
                while (abs(fxk - fxkp) > eps()) & (fxkp > fxk - qls_opt.gamma * alpha * alpha)
                    # if fxkp < fxk - qls_opt.gamma * alpha * d_flow2
                    alpha = alpha * qls_opt.delta
                    fxkp = lower_obj_p(flow + alpha * d_flow, td, path_ij, ep)
                end
            end

            # flow update and check if there are 0 flow path (which should be deleted)
            tmp.flows[ind_i] += alpha
            tmp.flows[ind_j] -= alpha
            # tmp.flows[ind_j] = max(tmp.flows[ind_j], 0)
            @assert minimum(tmp.flows[ind_j]) >= 0

            if tmp.flows[ind_i] == 0
                tmp.npath -= 1
                deleteat!(tmp.paths, ind_i)
                deleteat!(tmp.flows, ind_i)
                n_path_add -= 1
            elseif tmp.flows[ind_j] == 0
                tmp.npath -= 1
                deleteat!(tmp.paths, ind_j)
                deleteat!(tmp.flows, ind_j)
                n_path_del += 1
            end
            m_od[r][ss] = tmp

            flow += alpha * d_flow
            flow[path_ij[findall(x -> (x<0), flow[path_ij])]] .= 0
        end

    end

    return flow 
end

# optimize Phi w.r.t v
function ISMO_F!(m_od, Flow::Vector{Float64}, Flow0, td::TA_Data, ep, graph, link_dic, qls_opt, bi_para, require_generation = true)
    if require_generation
        local state::Graphs.DijkstraState{Float64, Int}
        # local state::Graphs.BellmanFordState{Float64, Int}
    end
    local flow_cost::Vector{Float64}
    n_path_add = 0
    n_path_del = 0
    
    for r in 1:size(td.travel_demand)[1]
        # for each origin node r, find shortest paths to all destination nodes
        flow_cost = dPhi(td, Flow, Flow0, ep, bi_para)

        if require_generation
            # state = TA_bf_shortest_paths(graph, flow_cost, r, td.init_node, td.term_node, td.first_thru_node)
            state = TA_dijkstra_shortest_paths(graph, flow_cost, r, td.init_node, td.term_node, td.first_thru_node)
        end

        # for s in 1:size(td.travel_demand)[2]
        for ss in 1:length(m_od[r])
            s = m_od[r][ss].D 
            # print(" s = ", s)
            # pass if no demand
            if td.travel_demand[r, s] == 0
                continue
            end
            
            # add path if needed 
            @views tmp = m_od[r][ss]

            if require_generation
                new_path = get_path(state, r, s, link_dic)
                if new_path == []
                    continue
                end
                if !(new_path in tmp.paths)
                    tmp.npath += 1
                    push!(tmp.paths, get_path(state, r, s, link_dic))
                    push!(tmp.flows, 0)
                    n_path_add += 1
                end
            end

            # pass if there is only 1 path for the OD pair
            if tmp.npath == 1; continue; end
            @assert tmp.npath > 0
            # calculate the direction
            ind_i, ind_j = 1, 1
            # cost_p = sum(flow_cost[tmp.paths[1]])
            cost_p = dPhi_p(td, Flow, Flow0, ep, tmp.paths[1], bi_para)
            cost_i = cost_p 
            cost_j = cost_p
            for k in 2:tmp.npath
                # cost_p = sum(flow_cost[tmp.paths[k]])
                cost_p = dPhi_p(td, Flow, Flow0, ep, tmp.paths[k], bi_para)
                if cost_p < cost_i 
                    ind_i = k 
                    cost_i = cost_p
                end
                if cost_p > cost_j 
                    ind_j = k 
                    cost_j = cost_p 
                end
            end
            if ind_i == ind_j; continue; end

            d_flow = zeros(td.number_of_links)
            d_flow[tmp.paths[ind_i]] .+= 1
            d_flow[tmp.paths[ind_j]] .-= 1
            @views beta_r = tmp.flows[ind_j]
            
            path_ij = setdiff(union(tmp.paths[ind_i], tmp.paths[ind_j]), intersect(tmp.paths[ind_i], tmp.paths[ind_j]))
            @assert length(path_ij) > 0
            # d_flow2 = length(path_ij)
            
            # QLS
            fxk = Phi_p(Flow, Flow0, td, ep, path_ij, bi_para)
            # alpha = beta_r
            alpha = min(beta_r, (cost_j - cost_i)/Phi_dHd(td, Flow, ep, path_ij, bi_para))
            # print("cal alpha; ")
            if alpha > eps()
                fxkp = Phi_p(Flow + alpha * d_flow, Flow0, td, ep, path_ij, bi_para)
                while (alpha > 2*eps()) & (abs(fxk - fxkp) > 2*eps()) & (fxkp > fxk - qls_opt.gamma * alpha * alpha)
                    alpha = alpha * qls_opt.delta
                    # @printf("%.4e\n", alpha)
                    fxkp = Phi_p(Flow + alpha * d_flow, Flow0, td, ep, path_ij, bi_para)
                    # println(beta_r, " ", (cost_j - cost_i)/Phi_dHd(td, path_ij, bi_para)," ", alpha)
                end
            end

            tmp.flows[ind_i] += alpha
            tmp.flows[ind_j] -= alpha
            @assert minimum(tmp.flows[ind_j]) >= 0

            if tmp.flows[ind_i] == 0
                tmp.npath -= 1
                deleteat!(tmp.paths, ind_i)
                deleteat!(tmp.flows, ind_i)
                n_path_add -= 1
            elseif tmp.flows[ind_j] == 0
                tmp.npath -= 1
                deleteat!(tmp.paths, ind_j)
                deleteat!(tmp.flows, ind_j)
                n_path_del += 1
            end
            m_od[r][ss] = tmp

            Flow += alpha * d_flow
            Flow[path_ij[findall(x -> (x<0), Flow[path_ij])]] .= 0
        end
    end

    return Flow
end

function enhance3!(ep::Enhancement_plan, y0, Flow::Vector{Float64}, flow::Vector{Float64}, td::TA_Data, bi_para, yU=Inf*ones(lastindex(y0)))
    for i in 1:lastindex(y0)
        enhance_ya!(ep, i, y0[i], Flow, flow, td, bi_para, yU[i])
    end
end

function enhance_ya!(ep::Enhancement_plan, ind, y0, Flow::Vector{Float64}, flow::Vector{Float64}, td::TA_Data, bi_para, yU = Inf)
    iind = ep.inds[ind]
    y1 = ep.y[ind]
    if td.power[iind] == 0
        ep.y[ind] = 0
    else
        dvk = -td.power[iind]/(td.power[iind] + 1) * td.free_flow_time[iind] * td.b[iind] * (flow[iind] / (td.capacity[iind] + y0))^(td.power[iind] + 1)

        for _ in 1:10
            tmp = td.free_flow_time[iind] * td.b[iind] * (Flow[iind] / (td.capacity[iind] + y1)) ^ (td.power[iind] + 1)
                
            dy = - (td.power[iind] + td.power[iind]/(td.power[iind] + 1) * bi_para.rho) * tmp + 2 * ep.d[ind] * y1 - bi_para.rho * dvk + bi_para.beta * (y1 - y0)

            if abs(dy) < 1e-5; break; end

            ddy = td.power[iind] * (td.power[iind] + 1 + bi_para.rho) * tmp / (td.capacity[iind] + y1) + 2 * ep.d[ind] + bi_para.beta 

            y1 = y1 - dy / ddy 
            if y1 < 0;  y1 = 0; break; end 
            if y1 > yU; y1 = yU; break; end
        end 
        ep.y[ind] = y1
    end
end

function enhance_S!(ep::Enhancement_plan, y0, Flow::Vector{Float64}, td::TA_Data, yU)
    for ind in 1:lastindex(y0)
        iind = ep.inds[ind]
        y1 = ep.y[ind]
        if td.power[iind] == 0
            ep.y[ind] = 0
        else
            for _ in 1:10
                tmp = td.free_flow_time[iind] * td.b[iind] * (Flow[iind] / (td.capacity[iind] + y1)) ^ (td.power[iind] + 1)
                    
                dy = - td.power[iind] * tmp + 2 * ep.d[ind] * y1 

                if abs(dy) < 1e-5; break; end

                ddy = td.power[iind] * (td.power[iind] + 1) * tmp / (td.capacity[iind] + y1) + 2 * ep.d[ind] 

                y1 = y1 - dy / ddy 
                if y1 < 0;  y1 = 0; break; end 
                if y1 > yU[ind]; y1 = yU[ind]; break; end
            end 
            ep.y[ind] = y1
        end
    end
end


# used in AL method
function AL_tg(flow, flow2, ep, td, mu, rho)
    @views inds = ep.inds
    @views d = ep.d
    y = zeros(td.number_of_links)
    y[inds] = ep.y
    
    tk = zeros(td.number_of_links)
    gk = zeros(lastindex(inds))

    h = lower_obj(td, flow, ep) - lower_obj(td, flow2, ep)
    coe1 = mu + rho * h
    for i in 1:td.number_of_links
        tk[i] = td.free_flow_time[i] * (1 + coe1) + td.free_flow_time[i] * td.b[i] * (1 + coe1 + td.power[i]) * (flow[i] / (td.capacity[i] + y[i])) ^ td.power[i]
        if tk[i] < 0; tk[i] = 0; end 
    end
    for i in 1:lastindex(inds) 
        ii = inds[i]
        tmp1 = (flow[ii] / (td.capacity[ii] + y[ii])) ^ (td.power[ii] + 1)
        tmp2 = (flow2[ii] / (td.capacity[ii] + y[ii])) ^ (td.power[ii] + 1)
        gk[i] = 2 * d[i] * y[ii] - td.power[ii] * td.free_flow_time[ii] * td.b[ii] * ((1 + coe1 / (td.power[ii] + 1)) * tmp1 - coe1 / (td.power[ii] + 1) * tmp2)
    end
    return tk, gk, h
end

function AL_INNER(flow, ep, yL, yU, td, al_para, graph, link_dic)
    mu = al_para.mu
    rho = al_para.rho
    EPS_IN = al_para.EPS_IN

    k = 0
    while true
        global h2 = 0
        k = k+1 
        # 1 DUE 
        flow2 = deepcopy(flow) 
        for i in 1:3
            flow2 = ta_fw(td, flow2, ep, graph, link_dic)
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
        # Zk = abs(tk'*d1 + gk'*d2)
        Zk = abs(tk'*d1 + gk'*d2) / (upper_obj(td, flow, ep) + mu*h2 + .5*rho*h2^2)
        # println(Zk)
        if Zk < EPS_IN; break; end 

        # 6 line search
        alpha = 1 / k

        # 7 update 
        flow = flow + alpha * d1 
        ep.y = ep.y + alpha * d2
        
        # @printf("upper_obj: %8.6f (aec = %5.2e), test = %5.2e\n", upper_obj(td, flow, dd, y), aec(flow, td, graph, link_dic, y), test)

    end 
    # @printf("     lower obj is %9.7f  new obj is %9.7f\n", lower_obj(td, flow, y), lower_obj(td, flow2, y))
    return flow, ep, h2
end
