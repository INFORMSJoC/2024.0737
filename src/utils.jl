# ==============================================================================
# This file contains a series of functions used for the CNDP solvers. 
# It includes functions for:
# - Initializing all-or-nothing assignments (all_or_nothing and all_or_nothing!)
# - Objective functions for optimization (ta_objective, lower_obj, upper_obj)
# - System optimal and link optimal flow calculations (dsys, dsys_p, sys_p)
# - Metrics for flow (aec, rgap, sysgap)
# - Methods for calculating xi for DC optimization and error metrics (cal_xi_dc, err_vy)
# - Functions for specific optimization methods (HJ method, EDO method)
# - Evaluation functions for testing and validating results (eval_ep_ISMO, Eval_EP_SiouxFall)
#
# Author: Haian Yin
# Date: 2022-09-06
# ==============================================================================

# Initialization
function all_or_nothing!(m_od, flow::Vector{Float64}, td::TA_Data, graph, link_dic)
    # used by ISMO
    local state::Graphs.DijkstraState{Float64, Int}

    for r in 1:size(td.travel_demand)[1]
        # for each origin node r, find shortest paths to all destination nodes
        state = TA_dijkstra_shortest_paths(graph, td.free_flow_time, r, td.init_node, td.term_node, td.first_thru_node)

        for s in 1:size(td.travel_demand)[2]
            # for each destination node s, find the shortest-path vector
            # add_demand_vector!(v, td.travel_demand[r, s], state, r, s, link_dic)
            if td.travel_demand[r, s] == 0
                continue
            end

            tmp = OD_sol(r, s, 1, [get_path(state, r, s, link_dic)], [td.travel_demand[r, s]])
            push!(m_od[r], tmp)
            # m_od[r,s] = tmp
            flow[tmp.paths[1]] .+= td.travel_demand[r, s]
        end
    end
end

function all_or_nothing(travel_time::Vector{Float64}, td::TA_Data, graph, link_dic)
    # used by frank-wolfe
    local state::Graphs.DijkstraState{Float64, Int}
    flow = zeros(td.number_of_links)

    for r in 1:size(td.travel_demand)[1]
        # for each origin node r, find shortest paths to all destination nodes
        state = TA_dijkstra_shortest_paths(graph, travel_time, r, td.init_node, td.term_node, td.first_thru_node)

        for s in 1:size(td.travel_demand)[2]
            # for each destination node s, find the shortest-path vector
            if td.travel_demand[r, s] == 0; continue; end
            add_demand_vector!(flow, td.travel_demand[r, s], state, r, s, link_dic)
        end
    end

    return flow 
end

# objectives
function ta_objective(x::Vector{Float64}, td::TA_Data)
    # value = free_flow_time .* ( x + b.* ( x.^(power+1)) ./ (capacity.^power) ./ (power+1))
    # return sum(value)

    sum = 0.0
    for i in 1:length(x)
        sum += td.free_flow_time[i] * ( x[i] + td.b[i]* ( x[i]^(td.power[i]+1)) / (td.capacity[i]^td.power[i]) / (td.power[i]+1))
        sum += td.toll_factor * td.toll[i] + td.distance_factor * td.link_length[i]
    end
    return sum
end

function dtdv(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    """
    calculate dt/dv
    """
    y = spzeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = zeros(td.number_of_links)
    for ind in 1:td.number_of_links
        if td.power[ind] > 0
            tmp[ind] = td.power[ind] * td.b[ind] * td.free_flow_time[ind] * flow[ind]^(td.power[ind] - 1) * (1 / (td.capacity[ind] + y[ind]))^(td.power[ind])
        end
    end
    tmp = tmp[flow.>0]
    return spdiagm(0 => tmp)
end

function dtdy(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    """
    dtdy = power * tau *0.15 * v^power * (1/(c+y))^(power+1)
    tau: td.free_flow_time
    v: flow
    c: td.capacity
    y: ep.y
    """
    y = spzeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = zeros(td.number_of_links, size(ep.inds)[1])
    for i in 1:td.number_of_links
        if (i in ep.inds) == true
            if td.power[i] > 0
                tmp[i, findall(x -> x == i, ep.inds)[1]] = -td.power[i] * td.b[i] * td.free_flow_time[i] * flow[i]^(td.power[i]) / (td.capacity[i] + y[i])^(td.power[i] + 1)
            end
        end
    end
    tmp = tmp[flow.>0, :]
    return tmp
end

function dphidy(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    """
    calculate dphi/dy, phi is F(y,S(y)), calculate the derivative of phi with respect to y
    """
    y = spzeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = zeros(td.number_of_links)
    for i in 1:td.number_of_links
        if td.power[i] > 0
            tmp[i] = -td.power[i] * td.b[i] * td.free_flow_time[i] * flow[i]^(td.power[i] + 1) / (td.capacity[i] + y[i])^(td.power[i] + 1)
        end
    end

    return tmp[ep.inds]
end

function get_gradients(dvy::Matrix{Float64}, td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan, update_method)
    """
    a function to calculate the gradient with chain rule 
        input:
            dvy: the gradient of the upper level objective function with respect to the flow
            td: the traffic data
            flow: the flow vector
            ep: the enhancement plan
            update_method: the method to update the upper level variable
        output:
            gradient: the gradient of the upper level objective function with respect to the enhancement plan
    """
    tmp1 = bpr(td, flow, ep)
    dtv = dtdv(td, flow, ep)
    if update_method == 0
        gradient = 2 * ep.d .* ep.y + dphidy(td, flow, ep)
    else
        gradient = dphidy(td, flow, ep)
    end
    tmp1 = tmp1[flow.>0, :]
    flow = flow[flow.>0, :]
    for i in 1:size(dvy)[1]
        gradient += (tmp1[i] + (dtv[i, i] * flow[i])) * dvy[i, :]
    end
    return gradient
end


function lower_obj(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = 0.0
    for i in 1:td.number_of_links
        if td.power[i] == 0
            tmp += td.free_flow_time[i] * (1 + td.b[i]) * flow[i]
        else 
            # tmp += td.free_flow_time[i] * flow[i] * (1 + td.b[i] .* (flow[i]/(td.capacity[i] + y[i])) ^ (td.power[i]) / (td.power[i]+1))
            tmp += td.free_flow_time[i] * (flow[i] + td.b[i] .* (flow[i]^(td.power[i]+1)/(td.capacity[i] + y[i])^td.power[i]) / (td.power[i]+1))
        end
        # tmp += ts.free_flow_time[i] * (ts.flow[i] + ts.b[i] .* ts.flow[i] ^ (ts.power[i] + 1) /ts.capacity_p[i] ^ (ts.power[i]) / (ts.power[i]+1))
    end
    return tmp
end

function upper_obj(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    # return F 
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = 0.0
    for i in 1:td.number_of_links
        if flow[i] > 0
            tmp += td.free_flow_time[i] * flow[i] * (1 + td.b[i] .* (flow[i]/(td.capacity[i] + y[i])) ^ (td.power[i]))
        end
    end
    return tmp + sum(ep.d .* ep.y .* ep.y)
end

function upper_objs(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    # return two terms in F
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = 0.0
    for i in 1:td.number_of_links
        if flow[i] > 0
            tmp += td.free_flow_time[i] * flow[i] * (1 + td.b[i] .* (flow[i]/(td.capacity[i] + y[i])) ^ (td.power[i]))
        end
    end
    return tmp, sum(ep.d .* ep.y .* ep.y)
end

# used for finding system optimal
function dsys(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = zeros(td.number_of_links)
    for i in 1:td.number_of_links
        flow[i] = max(flow[i], 0) 
        tmp[i] = td.free_flow_time[i] * (1 + (td.power[i] + 1) * td.b[i] * (flow[i] / (td.capacity[i] + y[i])) ^ td.power[i])
    end
    return tmp
end

function dsys_p(td::TA_Data, path, flow::Vector{Float64}, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = 0
    for i in path
        tmp += td.free_flow_time[i] * (1 + (td.power[i] + 1) * td.b[i] * (flow[i] / (td.capacity[i] + y[i])) ^ td.power[i])
    end 
    return tmp
end

function sys_p(flow::Vector{Float64}, td::TA_Data, path, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = 0
    for i in path
        flow[i] = max(flow[i], 0)
        if td.power[i]>0
            tmp += td.free_flow_time[i] * flow[i] * (1 + td.b[i] * (flow[i]/(td.capacity[i]+y[i])) ^ (td.power[i]))
        else 
            tmp += td.free_flow_time[i] * flow[i] * (1 + td.b[i]) 
        end
    end
    return tmp
end

# used for lower update
function bpr(td::TA_Data, flow::Vector{Float64}, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = zeros(td.number_of_links)
    for i in 1:td.number_of_links
        # flow[i] = max(flow[i], 0) 
        if td.power[i] == 0
            tmp[i] = td.free_flow_time[i] * (1.0 + td.b[i])
        else
            tmp[i] = td.free_flow_time[i] * (1.0 + td.b[i] * (flow[i] / (td.capacity[i] + y[i])) ^ td.power[i])
        end
    end
    return tmp
end

function bpr_p(td::TA_Data, path, flow::Vector{Float64}, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = 0
    for i in path
        tmp += td.free_flow_time[i] * (1 + td.b[i] * (flow[i] / (td.capacity[i] + y[i])) ^ td.power[i])
    end 
    return tmp
end

function lower_obj_p(flow::Vector{Float64}, td::TA_Data, path, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    lower_obj = 0
    for i in path
        flow[i] = max(flow[i], 0)
        if td.power[i]>0
            lower_obj += td.free_flow_time[i] * flow[i] * (1 + td.b[i] * (flow[i]/(td.capacity[i]+y[i])) ^ (td.power[i]) / (td.power[i]+1))
        else 
            lower_obj += td.free_flow_time[i] * flow[i] * (1 + td.b[i]) 
        end
    end
    return lower_obj
end

function lower_obj_p_dHd(flow::Vector{Float64}, td::TA_Data, path, ep::Enhancement_plan)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    lower_obj_p_dHd = 0
    for ind in path
        if td.power[ind] <= 1
            continue
        else 
            lower_obj_p_dHd += td.power[ind] * td.free_flow_time[ind] * td.b[ind] * (flow[ind] / (td.capacity[ind] + y[ind])) ^ (td.power[ind]-1) / (td.capacity[ind] + y[ind])
        end
    end
    return lower_obj_p_dHd
end

# used for upper update
function dPhi(td::TA_Data, Flow, flow0::Vector{Float64}, ep::Enhancement_plan, bi_para)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    tmp = zeros(td.number_of_links)
    for ind in 1:td.number_of_links
        tmp[ind] = (1+bi_para.rho) * td.free_flow_time[ind] + (td.power[ind] + 1 + bi_para.rho) * td.free_flow_time[ind] * td.b[ind] * (Flow[ind] / (td.capacity[ind] + y[ind])) ^ td.power[ind] + bi_para.beta * (Flow[ind] - flow0[ind])
    end
    return tmp
end

function dPhi_p(td::TA_Data, Flow::Vector{Float64}, flow0::Vector{Float64}, ep::Enhancement_plan, path, bi_para)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    res = 0
    for ind in path
        res += (1+bi_para.rho) * td.free_flow_time[ind] + ((td.power[ind]+1) + bi_para.rho) * td.free_flow_time[ind] * td.b[ind] * (Flow[ind] / (td.capacity[ind] + y[ind])) ^ td.power[ind] + bi_para.beta * (Flow[ind] - flow0[ind])
    end 
    return res 
end

function Phi_p(flow::Vector{Float64}, flow0::Vector{Float64}, td::TA_Data, ep::Enhancement_plan, path, bi_para)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    res = 0
    for ind in path
        if flow[ind]>0
            res += td.free_flow_time[ind] * flow[ind] * (1+bi_para.rho + (1+bi_para.rho/(td.power[ind]+1)) * td.b[ind] * (flow[ind] / (td.capacity[ind] + y[ind])) ^ td.power[ind]) + bi_para.beta/2 * (flow[ind] - flow0[ind])^2
        end
    end
    return res
end

function Phi_dHd(td::TA_Data, Flow::Vector{Float64}, ep::Enhancement_plan, path, bi_para)
    y = zeros(td.number_of_links)
    y[ep.inds] = ep.y
    res = 0
    for ind in path
        if td.power[ind] <= 1
            res += bi_para.beta
        else 
            res += td.power[ind] * (td.power[ind]+1+bi_para.rho) * td.free_flow_time[ind] * td.b[ind] * (Flow[ind] / (td.capacity[ind] + y[ind]))^(td.power[ind]-1) / (td.capacity[ind] + y[ind]) + bi_para.beta
        end
    end
    return res
end

# calculate aec
function aec(flow::Vector{Float64}, td::TA_Data, graph, link_dic, ep::Enhancement_plan)
    local state::Graphs.DijkstraState{Float64, Int}
    travel_time = bpr(td, flow, ep)
    flow0 = all_or_nothing(travel_time, td,  graph, link_dic)
    return (flow' * travel_time - flow0' * travel_time) / (sum(td.travel_demand))
end

# calculate rgap
function rgap(flow::Vector{Float64}, td::TA_Data, graph, link_dic, ep::Enhancement_plan)
    local state::Graphs.DijkstraState{Float64, Int}
    travel_time = bpr(td, flow, ep)
    flow0 = all_or_nothing(travel_time, td,  graph, link_dic)
    tmp = flow' * travel_time
    return (tmp - flow0' * travel_time) / (tmp)
end

function sysgap(flow::Vector{Float64}, td::TA_Data, graph, link_dic, ep::Enhancement_plan)
    local state::Graphs.DijkstraState{Float64, Int}
    travel_time = dsys(td, flow, ep)
    flow0 = all_or_nothing(travel_time, td,  graph, link_dic)
    tmp = flow' * travel_time
    return (tmp - flow0' * travel_time) / (tmp)
end

# cal xi for dc method 
function cal_xi_dc(td, Flow, Ep, flow, ep)
    @views inds = ep.inds
    vk = lower_obj(td, flow, ep)
    dvk = - td.power[inds] ./ (td.power[inds] .+ 1) .* td.free_flow_time[inds] .* td.b[inds] .* (flow[inds] ./ (td.capacity[inds] .+ ep.y)).^(td.power[inds] .+ 1)
    f = lower_obj(td, Flow, Ep)
    # xi = (f - vk - dvk' * (Ep.y - ep.y)) / f
    xi = (f - vk - dvk' * (Ep.y - ep.y))
    return xi
end

function cal_xi_dc(td, Flow, Ep, flow, ep, m_od_f, graph, link_dic, qls_opt)
    @views inds = ep.inds
    vk = lower_obj(td, flow, ep)
    dvk = - td.power[inds] ./ (td.power[inds] .+ 1) .* td.free_flow_time[inds] .* td.b[inds] .* (flow[inds] ./ (td.capacity[inds] .+ ep.y)).^(td.power[inds] .+ 1)
    f = lower_obj(td, Flow, Ep)
    xi = (f - vk - dvk' * (Ep.y - ep.y)) / f
    # xi = (f - vk - dvk' * (Ep.y - ep.y))
    xi_iter = 0
    while xi < 0
        @printf("f is %.4e vk is %.4e xi is %.4e aec is %.4e and aec is %.4e\n", f, vk, xi, rgap(flow, td, graph, link_dic, ep), rgap(Flow, td, graph, link_dic, Ep))
        xi_iter += 1
        if xi_iter > 100; xi = 0; break; end
        flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
        # flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt, false)
        f = lower_obj(td, flow, Ep)
        vk = lower_obj(td, flow, ep)
        dvk = - td.power[inds] ./ (td.power[inds] .+ 1) .* td.free_flow_time[inds] .* td.b[inds] .* (flow[inds] ./ (td.capacity[inds] .+ ep.y)).^(td.power[inds] .+ 1)
        xi = (f - vk - dvk' * (Ep.y - ep.y)) / f
    end
    # CSV.write("flow.csv", DataFrame(flow; auto=true))
    # CSV.write("y.csv", DataFrame(ep.y; auto=true))
    return xi
end

function err_vy(flow, y, flow0, y0)
    tmp1 = sum((flow - flow0).^2)
    tmp2 = sum((y - y0).^2)
    # return sqrt((tmp1 + tmp2) / (sum(flow.^2) + sum(y.^2)))
    # return sqrt(tmp1 + tmp2) 
    return sqrt(tmp2)
    # return sqrt(tmp2/sum(y.^2))
end

# for hj method 
function HJ_Fz(td, flow, ind, d, z)
    return td.free_flow_time[ind] * flow[ind] * (1 + td.b[ind] .* (flow[ind]/(td.capacity[ind] + z)) ^ (td.power[ind])) + d * z^2
end

function HJ_Fy(td, flow, ep)
    @views inds = ep.inds 
    @views d = ep.d
    @views y = ep.y 
    yy = zeros(td.number_of_links)
    yy[inds] = y
    tmp = 0.0
    for i in 1:td.number_of_links
        if flow[i] > 0
            tmp += td.free_flow_time[i] * flow[i] * (1 + td.b[i] .* (flow[i]/(td.capacity[i] + yy[i])) ^ (td.power[i]))
        end
    end

    return tmp + sum(d .* y .* y)
end

# for EDO method
function EDO_dZa(td, ind, d_a, flow_a, y_a)
    return 2 * d_a * y_a - td.power[ind] * td.free_flow_time[ind] * td.b[ind] * (flow_a / (td.capacity[ind] + y_a)) ^ (td.power[ind] + 1)
end

# eval ep 
function eval_ep_ISMO(td, ep)
    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    m_od = Array{Any}(nothing, size(td.travel_demand)[1], size(td.travel_demand)[2])
    flow = zeros(td.number_of_links)

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)

    all_or_nothing!(m_od, flow, td, graph, link_dic);
    for i in 1:100
        if i % 5 == 1
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        else
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
        end
    end
    result = upper_obj(td, flow, ep)
    @printf("Upper Obj is %8.6f (rgap = %8.2e)\n", result, rgap(flow, td, graph, link_dic, ep))
end

# for evaluation 
function Eval_EP_SiouxFall(y; jam = 1)
    td, flow, ep = load_SiouxFalls(jam)
    inds = [16, 17, 19, 20, 25, 26, 29, 39, 48, 74]
    d =  0.001*[26., 40., 26., 40., 25., 25., 48., 34., 48., 34.]
    ep = Enhancement_plan(y, d, inds)
    flow = zeros(td.number_of_links)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end

    # Initialization
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    while rgap(flow, td, graph, link_dic, ep) > 1e-10
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
    end

    @printf("Upper Obj is %8.6f (rgap = %8.4e)\n", upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))

end
