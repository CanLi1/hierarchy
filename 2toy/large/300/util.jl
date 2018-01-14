mutable struct node 
    mult_p1
    mult_p2
    mult_x1
    mult_x2
    g
    leftnode
    rightnode
    vars_lb 
    vars_ub 
    p1_opt
    p2_opt
    x1_opt
    x2_opt
    v_node
    V_node
    fathom_type
    master_time 
    sub_time
    ubsub_time
end
node() = node([],[],[],[],[],[],[],[0.0 0.0],[4.0 2.0],0.0,0.0,0.0,0.0,0.0,0.0, :notfathomed, [], [], [])

function psolve(d_input)
    m = generate_sub(p1bar=d_input[:p1], p2bar=d_input[:p2], x1bar=d_input[:x1], x2bar=d_input[:x2], demand=d_input[:demand], prob=d_input[:prob], lb = d_input[:vars_lb], ub=d_input[:vars_ub])
    status = solve(m)
    if status != :Optimal
        m = generate_mosek_sub(p1bar=d_input[:p1], p2bar=d_input[:p2], x1bar=d_input[:x1], x2bar=d_input[:x2], demand=d_input[:demand], prob=d_input[:prob], lb = d_input[:vars_lb], ub=d_input[:vars_ub])
        status = solve(m)
    end
    if status != :Optimal
        m = generate_knitro_sub(p1bar=d_input[:p1], p2bar=d_input[:p2], x1bar=d_input[:x1], x2bar=d_input[:x2], demand=d_input[:demand], prob=d_input[:prob], lb = d_input[:vars_lb], ub=d_input[:vars_ub])
        status = solve(m)
    end    
    if status != :Optimal
        error("sub not solved to optimality ")
    end

    d_output = Dict()
    d_output[:obj] = getobjectivevalue(m)
    d_output[:dual_p1] = getdual(getindex(m, :t1))
    d_output[:dual_p2] = getdual(getindex(m, :t2))
    d_output[:dual_x1] = getdual(getindex(m, :t3))
    d_output[:dual_x2] = getdual(getindex(m, :t4))
    d_output[:g] = d_output[:obj] - d_input[:p1] * d_output[:dual_p1]- d_input[:p2] * d_output[:dual_p2] - d_input[:x1] * d_output[:dual_x1]- d_input[:x2] * d_output[:dual_x2]

    return d_output
end

function psolve_ub(d_input)
    m = generate_ubsub(p1bar=d_input[:p1], p2bar=d_input[:p2], x1bar=d_input[:x1], x2bar=d_input[:x2], demand=d_input[:demand], prob=d_input[:prob])
    status = solve(m)
    d_output = Dict()   
    d_output[:obj] = getobjectivevalue(m)
    if status != :Optimal
        error("ub subproblem not solved to optimality")
    end
    return d_output
end

function solve_node(n::node)
    a = now()
    b=now()
    master_time = a-b
    sub_time = a-b 
    ubsub_time = a-b 
    ub = 1e3
    lb = -1e3
    ub_record = []
    lb_record = []
    mult_p1 = n.mult_p1
    mult_p2 = n.mult_p2
    mult_x1 = n.mult_x1
    mult_x2 = n.mult_x2
    g = n.g

    if length(mult_p1) == 0
        m = generate_master()
    else
        m = generate_master(mult_p1=mult_p1, mult_p2=mult_p2,mult_x1=mult_x1, mult_x2=mult_x2, g=g, lb=n.vars_lb, ub=n.vars_ub, iter=1:length(mult_p1[1,:]))
    end

    while ub >= lb +1e-3
        #solve master problem 
        a = now()
        solve(m)
        b=now()
        master_time = master_time + b - a
        global p1 = getvalue(getindex(m, :p1))
        global p2 = getvalue(getindex(m, :p2))
        global x1 = getvalue(getindex(m, :x1))
        global x2 = getvalue(getindex(m, :x2))

        #change the first stage decisions and solve subproblem 
        temp_mult_p1 = zeros(length(scenarios), 1)
        temp_mult_p2 = zeros(length(scenarios), 1)
        temp_mult_x1 = zeros(length(scenarios), 1)
        temp_mult_x2 = zeros(length(scenarios), 1)
        temp_mult_g = zeros(length(scenarios), 1)
        temp_sub_stat = []
        temp_ub_stat = []
        sub_problem= []
        for i in scenarios
            d_input = Dict()
            d_input[:p1] = p1
            d_input[:p2] = p2
            d_input[:x1] = x1 
            d_input[:x2] = x2 
            d_input[:demand] = d[i]
            d_input[:prob] = prob[i]
            d_input[:vars_lb] = n.vars_lb
            d_input[:vars_ub] = n.vars_ub
            push!(sub_problem, d_input)
        end
        a = now()
        solve_output = pmap(psolve, sub_problem)
        b = now()
        sub_time = sub_time + b - a
        for i in scenarios
            temp_mult_p1[i] = solve_output[i][:dual_p1]
            temp_mult_p2[i] = solve_output[i][:dual_p2]
            temp_mult_x1[i] = solve_output[i][:dual_x1]
            temp_mult_x2[i] = solve_output[i][:dual_x2]
            temp_mult_g[i] = solve_output[i][:g]
        end


        if length(mult_p1) == 0
            mult_p1 = temp_mult_p1
            mult_p2 = temp_mult_p2
            mult_x1 = temp_mult_x1
            mult_x2 = temp_mult_x2
            g = temp_mult_g
            
        else
            mult_p1 = [mult_p1 temp_mult_p1]
            mult_p2 = [mult_p2 temp_mult_p2]
            mult_x1 = [mult_x1 temp_mult_x1]
            mult_x2 = [mult_x2 temp_mult_x2]        
            g = [g temp_mult_g]
        end 


        #record objective
        temp_ub = 0.0
        for i in scenarios
            temp_ub = temp_ub + solve_output[i][:obj]
        end
        push!(ub_record, temp_ub)
        push!(lb_record, getobjectivevalue(m))
        if lb < lb_record[end]
            lb = lb_record[end]
        end
        if ub > ub_record[end]
            ub = ub_record[end]
        end

        #generate new master problem
        m = generate_master(mult_p1=mult_p1, mult_p2=mult_p2,mult_x1=mult_x1, mult_x2=mult_x2, g=g, lb=n.vars_lb, ub=n.vars_ub, iter=1:length(mult_p1[1,:]))
    end 

    #solve m again to obtain the optimal fist stage solution and lower bound
    solve(m)
    n.x1_opt = getvalue(getindex(m, :x1))
    n.x2_opt = getvalue(getindex(m, :x2))
    n.p1_opt = getvalue(getindex(m, :p1))
    n.p2_opt = getvalue(getindex(m, :p2))
    n.v_node = getobjectivevalue(m)

    #update multipliers
    n.mult_p1 = mult_p1
    n.mult_p2 = mult_p2
    n.mult_x1 = mult_x1
    n.mult_x2 = mult_x2
    n.g = g 

    #update upper bound 
    ubsub_problem= []
    for i in scenarios
        d_input = Dict()
        d_input[:p1] = n.p1_opt
        d_input[:p2] = n.p2_opt
        d_input[:x1] = n.x1_opt 
        d_input[:x2] = n.x2_opt 
        d_input[:demand] = d[i]
        d_input[:prob] = prob[i]
        push!(ubsub_problem, d_input)
    end
    a = now()
    solve_output = pmap(psolve_ub, ubsub_problem)
    b = now()
    ubsub_time = ubsub_time + b - a 
    V_node = 0.0
    for i in scenarios
        V_node = V_node + solve_output[i][:obj]
    end
    n.V_node = V_node

    n.ubsub_time = ubsub_time
    n.sub_time = sub_time
    n.master_time = master_time
end












