# Large 
MAX_TIME = 1000

## DC 
function DC_tntp(root, name, U, per, seed, dmean, dstd, dc_para, jam = 1, due_acc = 1e-2, verbose=false)
    td, flow, ep = load_tntp(name, per, seed, dmean, dstd, jam)
    yU = U * td.capacity[ep.inds]
    ep.y = yU / 2

    name = name * "_" * string(jam) * "_"
    open(root * "DC_" * name * string(per) * ".txt", "w") do io; end
    open(root * "DC_" * name * string(per) * "_y.txt", "w") do io; end

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    m_od_f = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od_f[i] = Array{Any}(nothing, 0)
    end

    # Initialization
    T = 0
    timer = time()
    all_or_nothing!(m_od_f, flow, td, graph, link_dic);
    nii = 0
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        nii += 1
        flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
        if nii > 50; break; end
    end

    m_od_F = deepcopy(m_od_f)
    Flow = deepcopy(flow)
    T += time() - timer

    open(root * "DC_" * name * string(per) * ".txt", "a") do io
        @printf(io, "%8.2f %18.2f %9.2e\n", T, upper_obj(td, Flow, ep), rgap(Flow, td, graph, link_dic, ep))
    end

    while true
        timer = time()
        ep0 = deepcopy(ep)
        Flow0 = deepcopy(Flow)
        
        # solving subproblem1
        for _ in 1:1
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
        end
        while rgap(flow, td, graph, link_dic, ep) > due_acc
            flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
        end

        # solving subproblem2
        for _ in 1:1
            enhance3!(ep, ep0.y, Flow, flow, td, dc_para, yU)
            Flow = ISMO_F!(m_od_F, Flow, Flow0, td, ep, graph, link_dic, qls_opt, dc_para)
            enhance3!(ep, ep0.y, Flow, flow, td, dc_para, yU)
        end

        xi = cal_xi_dc(td, Flow, ep, flow, ep0, m_od_f, graph, link_dic, qls_opt)
        err = err_vy(Flow, ep.y, Flow0, ep0.y)

        T += time() - timer
        open(root * "DC_" * name * string(per) * ".txt", "a") do io
            @printf(io, "%8.2f %18.2f %18.2f %9.2e %9.2e %9.2e %9.2e %6d\n", T, upper_obj(td, Flow, ep), upper_obj(td, flow, ep), rgap(Flow, td, graph, link_dic, ep), xi / lower_obj(td, Flow, ep), xi, err, dc_para.rho)
        end

        if xi / lower_obj(td, Flow, ep) < dc_para.EPS; break; end
        if max(dc_para.rho, 1/xi) < maximum(yU)/err
            # dc_para.rho += dc_para.tau
            dc_para.rho *= dc_para.tau
        end

        if T > MAX_TIME; break; end 

    end

    open(root * "DC_" * name * string(per) * "_y.txt", "a") do io
        for i in 1:lastindex(ep.y); @printf(io, "%6.4f ", ep.y[i]); end
        @printf(io, "\n")
    end
end

## AL 
function AL_tntp(root, name, U, per, seed, dmean, dstd, al_para, jam = 1, due_acc = 1e-2, verbose=false)
    td, flow, ep = load_tntp(name, per, seed, dmean, dstd, jam)
    name = name * "_" * string(jam) * "_"
    open(root * "AL_" * name * string(per) * ".txt", "w") do io; end
    open(root * "AL_" * name * string(per) * "_y.txt", "w") do io; end

    # AL 
    @views inds = ep.inds
    yL = zeros(lastindex(inds))
    yU = U * td.capacity[inds]
    ep.y = yU / 2

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
    all_or_nothing!(m_od, flow, td, graph, link_dic)
    # for i in 1:10
    nii = 0
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        nii += 1
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        if nii > 50; break; end
    end

    flow2 = deepcopy(flow)
    for _ in 1:1
        flow2 = ISMO_f!(m_od, flow2, td, ep, graph, link_dic, qls_opt)
    end
    h0 = lower_obj(td, flow, ep) - lower_obj(td, flow2, ep)
    h2 = 0
    T += time() - timer

    open(root * "AL_" * name * string(per) * ".txt", "a") do io
        @printf(io, "%8.2f %18.2f %9.2e \n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
    end

    mu = al_para.mu
    rho = al_para.rho
    EPS_IN = al_para.EPS_IN
    # F0 = upper_obj(td, flow, ep) + mu * h0 + .5 * rho * h0^2
    F0 = 0
    k = 0
    for i = 1:100
        # global k
        # 1 inner loop 
        mu = al_para.mu
        rho = al_para.rho
        EPS_IN = al_para.EPS_IN

        
        # while k < 10000
        while true
            timer = time()
            # global h2 = 0
            # h2 = 0
            # global k = k+1 
            k = k+1
            # 1 DUE 
            for i in 1:1
                flow2 = ISMO_f!(m_od, flow2, td, ep, graph, link_dic, qls_opt)
            end
        
            # 2 calculation marginal    
            # 3 calculation general cost 
            tk, gk, h2 = AL_tg(flow, flow2, ep, td, mu, rho)

            # @assert minimum(tk) >= 0 "the minimum is $(minimum(tk))"
            
            # 4 all or nothing for v / y 
            flow3 = all_or_nothing(tk, td, graph, link_dic)
            y3 = zeros(lastindex(ep.inds))
            for i in 1:lastindex(ep.inds)
                if gk[i] > 0; y3[i] = yL[i]
                else; y3[i] = yU[i]
                end
            end

            # 5 stopping test
            d1 = flow3 - flow 
            d2 = y3 - ep.y 
            # if F0 == 0; F0 = abs(tk' * d1 + gk' * d2); end
            # Zk = abs(tk'*d1 + gk'*d2) / F0
            Zk = abs(tk'*d1 + gk'*d2) / (upper_obj(td, flow, ep) + mu * h2 + .5 * rho * h2^2)
            # Zk = abs(tk' * d1 + gk'*d2) / abs(tk' * flow3 + gk' * y3)
            # F1 = upper_obj(td, flow, ep) + mu * h2 + .5 * rho * h2^2
            # Zk = abs(F1 - F0) / F0
            if Zk < EPS_IN
                T += time() - timer
                open(root * "AL_" * name * string(per) * ".txt", "a") do io
                    @printf(io, "%8.2f %18.2f %9.2e %9.2e %9.2e out   %9.2e %9.2e\n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep), h2 / lower_obj(td, flow, ep), Zk, mu, rho)
                end
                break
            end 

            # F0 = F1

            # 6 line search
            # alpha = al_para.alpha0 / k
            # alpha = 10 / (k + 100)
            alpha = al_para.alpha1 / (k + al_para.alpha2)

            # 7 update 
            flow = flow + alpha * d1 
            ep.y = ep.y + alpha * d2
            T += time() - timer
            open(root * "AL_" * name * string(per) * ".txt", "a") do io
                @printf(io, "%8.2f %18.2f (%9.2e %9.2e) %9.2e %9.2f inner (%d) %9.2f %9.2e\n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep), rgap(flow2, td, graph, link_dic, ep), h2/lower_obj(td, flow, ep), Zk, k, mu, rho)
            end

            if T > MAX_TIME; break; end 
            
        end 

        open(root * "AL_" * name * string(per) * ".txt", "a") do io
            @printf(io, "h0 = %9.2e h2 = %9.2e\n", h0, h2)
        end
        # 2 stopping test
        # println(h2 / lower_obj(td, flow, ep))
        if h2 / lower_obj(td, flow, ep) < al_para.EPS_AL
            break
        end

        if T > MAX_TIME; break; end 

        # 3 updating Lagrangian multiplie
        # 4 update penalty parameter
        if h2 < 0.25 * h0
            al_para.mu = al_para.mu + al_para.rho * h2
        else        
            al_para.rho = al_para.beta * al_para.rho 
            if al_para.rho > 1e15; break; end
        end
        h0 = h2

    end

    # @printf("print y: \n")
    open(root * "AL_" * name * string(per) * "_y.txt", "w") do io
        for i in 1:lastindex(ep.y); @printf(io, "%6.4f ", ep.y[i]); end
        @printf(io, "\n")
    end
end


## HJ 
function HJ_tntp(root, name, U, per, seed, dmean, dstd, para, jam = 1, due_acc = 1e-2, verbose=false)
    td, flow, ep = load_tntp(name, per, seed, dmean, dstd, jam)
    name = name * "_" * string(jam) * "_"
    open(root * "HJ_" * name * string(per) * ".txt", "w") do io; end
    open(root * "HJ_" * name * string(per) * "_y.txt", "w") do io; end
    @views inds = ep.inds
    @views d = ep.d
    # para = MutableNamedTuple(eps = 1e-1, alpha = 0.9, delta = delta)

    yU = U * td.capacity[inds]
    ep.y = yU / 2
    # print(yU)
    # return 

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    # m_od = Array{Any}(nothing, size(td.travel_demand)[1], size(td.travel_demand)[2])
    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
 
    # Initialization
    T = 0
    timer = time()
    all_or_nothing!(m_od, flow, td, graph, link_dic);
    nii = 0
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        nii += 1
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        if nii > 50; break; end
    end

    # HJ method
    iter = 0
    z = deepcopy(ep.y)
    # K = div(length(inds), 10)
    T += time() - timer
    open(root * "HJ_" * name * string(per) * ".txt", "a") do io
        @printf(io, "%8.2f %18.2f %9.2e\n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
    end

    while true
        timer = time()
        iter += 1
        m_od_z = deepcopy(m_od)
        flow_z = deepcopy(flow)
        ep_z = deepcopy(ep)
        T += time() - timer
        
        for i in 1:lastindex(inds)
            timer = time()
            tmp = HJ_Fz(td, flow_z, inds[i], d[i], z[i])
            if HJ_Fz(td, flow_z, inds[i], d[i], min(z[i] + para.delta, yU[i])) < tmp
                z[i] = min(z[i] + para.delta, yU[i])
                ep_z.y[i] = z[i]
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt)
            elseif HJ_Fz(td, flow_z, inds[i], d[i], max(z[i] - para.delta, 0)) < HJ_Fz(td, flow_z, inds[i], d[i], z[i])
                z[i] = max(z[i] - para.delta, 0)
                ep_z.y[i] = z[i]
                flow_z = ISMO_f!(m_od_z, flow_z, td, ep_z, graph, link_dic, qls_opt)
            end 
            T += time() - timer
            if T > MAX_TIME; break; end 
        end

        timer = time()
        nii = 0

        if HJ_Fy(td, flow_z, ep_z) < HJ_Fy(td, flow, ep)
            y0 = deepcopy(ep.y)
            ep.y = deepcopy(z)
            nii = 0

            # while rgap(flow, td, graph, link_dic, ep) > due_acc
                # nii += 1
            for _ in 1:3
                flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
                # if nii > 50; break; end
            end         
            z = y0 + para.alpha * (ep.y - y0)
        elseif para.delta <= para.eps
            break 
        else 
            para.delta *= 0.5
            z = deepcopy(ep.y)
        end

        T += time() - timer
        open(root * "HJ_" * name * string(per) * ".txt", "a") do io
            @printf(io, "%8.2f %18.2f %9.2e\n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
        end

        if T > MAX_TIME; break; end 
    end

    open(root * "HJ_" * name * string(per) * "_y.txt", "a") do io
        for i in 1:lastindex(ep.y); @printf(io, "%6.4f ", ep.y[i]); end
        @printf(io, "\n")
    end

end


## EDO
function EDO_tntp(root, name, U, per, seed, dmean, dstd, EPS_EDO, jam = 1, due_acc = 1e-2, verbose=false)
    td, flow, ep = load_tntp(name, per, seed, dmean, dstd, jam)
    name = name * "_" * string(jam) * "_"
    open(root * "EDO_" * name * string(per) * ".txt", "w") do io; end
    open(root * "EDO_" * name * string(per) * "_y.txt", "w") do io; end

    # EDO 
    @views inds = ep.inds
    yL = zeros(lastindex(inds))
    # yU = yU * ones(lastindex(inds))
    yU = U * td.capacity[inds]
    ep.y = yU / 2

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
    nii = 0
    while rgap(flow, td, graph, link_dic, ep) > due_acc
        nii += 1
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        if nii > 50; break; end
    end
    T += time() - timer
    open(root * "EDO_" * name * string(per) * ".txt", "a") do io
        @printf(io, "%8.2f %18.2f %9.2e\n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
    end

    # EDO 
    timer = time()
    flag = true
    @views inds = ep.inds
    @views d = ep.d
    yy = (yL + yU) / 2
    iter = 0
    T += time() - timer

    while flag 
        flag = false
        iter += 1
        for i in 1:lastindex(inds)
            timer = time()
            if ~flag 
                if yU[i] - yL[i] < EPS_EDO 
                    continue
                else
                    flag = true
                end
            end

            yy[i] = (yL[i] + yU[i]) / 2
            tmp = EDO_dZa(td, inds[i], d[i], flow[inds[i]], yy[i])

            if tmp < 0
                yL[i] = yy[i] 
            elseif tmp > 0 
                yU[i] = yy[i] 
            else 
                yL[i] = yy[i]
                yU[i] = yy[i]
            end
            T += time() - timer

            if T > MAX_TIME 
                open(root * "EDO_" * name * string(per) * ".txt", "a") do io
                    @printf(io, "%8.2f %18.2f %9.2e\n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
                end
                break 
            end 
        end

        timer = time()
        ep.y = yy

        # solving subproblem
        for _ in 1:5
            flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        end
        
        nii = 0

        T += time() - timer
        open(root * "EDO_" * name * string(per) * ".txt", "a") do io
            @printf(io, "%8.2f %18.2f %9.2e\n", T, upper_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
        end

        if T > MAX_TIME; break; end 
    end
    
    open(root * "EDO_" * name * string(per) * "_y.txt", "a") do io
        for i in 1:lastindex(ep.y); @printf(io, "%6.4f ", ep.y[i]); end
        @printf(io, "\n")
    end

end

## Post processing
function post_tntp(root, network_name, seed, dmean, dstd, eps, U = 1, jam=1, methods = ["DC", "AL", "EDO", "HJ"])

    td = load_ta_network(network_name)
    td.travel_demand *= jam
    name = network_name * "_" * string(jam) * "_"
    open(root * "post" * name * ".txt", "w") do io; end

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    m_od = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od[i] = Array{Any}(nothing, 0)
    end
    flow = zeros(td.number_of_links)

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)

    all_or_nothing!(m_od, flow, td, graph, link_dic);

    ep = Enhancement_plan([0.], [0.], [1])

    i = 0
    open(root * "post" * name * ".txt", "a") do io
        @printf(io, "\n\n %4s Origin:\n", method)
    end

    while (rgap(flow, td, graph, link_dic, ep) > eps) & (i<101)
        i += 1
        flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
        open(root * "post" * name * ".txt", "a") do io
            @printf(io, "%3d-th iteration, upper_obj: %18.6f, lower_obj: %18.6f, rgap: %10.4e\n", i, upper_obj(td, flow, ep), lower_obj(td, flow, ep), rgap(flow, td, graph, link_dic, ep))
        end
    end
    flow0 = deepcopy(flow)
    m_od0 = deepcopy(m_od)

    for per in [.1, .3, .7]
        open(root * "post" * name * ".txt", "a") do io
            @printf(io, "\n\n per = %3.1f\n", per)
        end
        rng = MersenneTwister(seed)
        inds = randsubseq(rng, 1:td.number_of_links, per)
        d = 1e-3 .* (dmean .+ dstd .*rand(rng, length(inds)))

        for method in methods
            open(root * "post" * name * ".txt", "a") do io
                @printf(io, "\n\n %4s METHOD:\n", method)
            end

            ydata = readdlm(root * "" * method * "_" * name * string(per) * ".txt")
            if ydata[end, 1] == "h0"
                time = ydata[end-1, 1]
            else
                time = ydata[end, 1]
            end
            ydata = readdlm(root * "" * method * "_" * name * string(per) * "_y.txt")
            y = ydata[end, :]
            open(root * "post" * name * ".txt", "a") do io
                @printf(io, "\n\n Time: %.2f\n", time)
            end

            ep = Enhancement_plan(y, d, inds)

            F0, f0 = upper_obj(td, flow, ep), lower_obj(td, flow, ep)
            open(root * "post" * name * ".txt", "a") do io
                @printf(io, "%3d-th iteration, upper_obj: %18.6f, lower_obj: %18.6f, rgap: %10.4e\n", 0, F0, f0, rgap(flow, td, graph, link_dic, ep))
            end
            flow = deepcopy(flow0)
            m_od = deepcopy(m_od0)


            for i in 1:100
                if i % 5 == 1
                    flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt)
                else
                    flow = ISMO_f!(m_od, flow, td, ep, graph, link_dic, qls_opt, false)
                end

                F1 = upper_obj(td, flow, ep)
                f1 = lower_obj(td, flow, ep)

                open(root * "post" * name * ".txt", "a") do io
                    @printf(io, "%3d-th iteration, upper_obj: %18.6f, lower_obj: %18.6f, rgap: %10.4e\n", i, F1, f0, rgap(flow, td, graph, link_dic, ep))
                end
                if abs(f1 - f0) / f1 < 1e-5
                    break
                else
                    F0, f0 = F1, f1
                end
            end

        end

    end

end 
