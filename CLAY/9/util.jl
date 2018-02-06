mutable struct node 
	mult_x
	mult_y
    g
    leftnode
    rightnode
    xlb
    xub
    ylb
    yub
    x_opt
    y_opt
    v_node
    V_node
    fathom_type
    master_time 
    sub_time
    ubsub_time
end
node() = node([],[],[],[],[], xlb, xub, ylb, yub, xlb, ylb, -1e8, 1e8, :notfathomed, [], [], [])

function psolve(m::JuMP.Model)
	d_output = Dict()
	d_output[:error] = false
	try
    	global status = solve(m)
    catch err 
    	d_output[:error] = true
    end

    if status != :Optimal
        d_output[:error] = true
    end

    if !d_output[:error]
	    d_output[:obj] = getobjectivevalue(m)
	    d_output[:dual_x] = getdual(getindex(m, :t1))
	    d_output[:dual_y] = getdual(getindex(m, :t2))
	    d_output[:model] = m
    end
    return d_output
end


# function psolve_ub(d_input)
#     m = generate_ubsub(xbar=d_input[:xbar], ybar=d_input[:ybar], Carea=d_input[:Carea], prob = d_input[:prob])
#     status = solve(m)
#     d_output = Dict()   
#     d_output[:obj] = getobjectivevalue(m)
#     if status != :Optimal
#         error("ub subproblem not solved to optimality")
#     end
# 	return d_output
# end


function solve_node(n::node)
    a = now()
    b=now()
    master_time = a-b
    sub_time = a-b 
    ubsub_time = a-b 
    ub = 1e8
    lb = -1e8
    ub_record = []
    lb_record = []
    mult_x = n.mult_x 
    mult_y = n.mult_y
    g = n.g

    if length(mult_x) == 0
        m = generate_master(xub=n.xub, xlb=n.xlb, ylb=n.ylb, yub=n.yub)
    else
        m = generate_master(mult_x=mult_x, mult_y=mult_y, g=g, xlb=n.xlb, xub=n.xub, ylb =n.ylb, yub=n.yub, iter=1:length(mult_x))
    end

    #generate sub_problems
    sub_problem= []
    for s in scenarios
    	push!(sub_problem, generate_dnf(Carea=Carea[:,s], prob=prob[s], xub=n.xub, xlb=n.xlb, ylb=n.ylb, yub=n.yub))
    end

    while ub >= lb * 1.0001
        #solve master problem 
        a = now()
        solve(m)
        b=now()
        master_time = master_time + b - a
        global x = getvalue(getindex(m, :x))
        global y = getvalue(getindex(m, :y))

        #change the first stage decisions and solve subproblem 
        temp_mult_x = []
        temp_mult_y = []
        temp_g = []
        
        for s in scenarios
        	push!(temp_mult_x, zeros(length(rectangles)))
        	push!(temp_mult_y, zeros(length(rectangles)))
        	push!(temp_g, 0.0)

        	#update first stage solutions
        	for i in rectangles
        		JuMP.setRHS(getindex(sub_problem[s], :t1)[i], x[i])
        		JuMP.setRHS(getindex(sub_problem[s], :t2)[i], y[i])
        	end
        end
            
       
        a = now()
        solve_output = pmap(psolve, sub_problem)
        b = now()
        sub_time = sub_time + b - a

        #resolve again using knitro if ipopt fails
        for s in scenarios
        	if solve_output[s][:error]
        		sub_problem[s] = generate_knitro_dnf(xbar=x, ybar=y, Carea=Carea[:,s], prob=prob[s], xub=n.xub, xlb=n.xlb, ylb=n.ylb, yub=n.yub)
        		solve(sub_problem[s])
			    solve_output[s][:obj] = getobjectivevalue(sub_problem[s])
			    solve_output[s][:dual_x] = getdual(getindex(sub_problem[s], :t1))
			    solve_output[s][:dual_y] = getdual(getindex(sub_problem[s], :t2))  
			    sub_problem[s] = generate_dnf(xbar=x, ybar=y, Carea=Carea[:,s], prob=prob[s], xub=n.xub, xlb=n.xlb, ylb=n.ylb, yub=n.yub)
			else
				sub_problem[s] = solve_output[s][:model]
			end
		end

        for s in scenarios
        	for i in rectangles
        		temp_mult_x[s][i] = solve_output[s][:dual_x][i]
        		temp_mult_y[s][i] = solve_output[s][:dual_y][i]
        	end
        	temp_g[s] = solve_output[s][:obj] - sum(solve_output[s][:dual_x][i] * x[i] + solve_output[s][:dual_y][i] * y[i] for i in rectangles)
        end

        #update multipliers
        push!(mult_x, temp_mult_x)
        push!(mult_y, temp_mult_y)
        push!(g, temp_g)

        #record objective
        temp_ub = getvalue(getindex(m, :firststagecost))
        for s in scenarios
            temp_ub = temp_ub + solve_output[s][:obj]
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
        m = generate_master(mult_x=mult_x, mult_y=mult_y, g=g, xlb=n.xlb, xub=n.xub, ylb =n.ylb, yub=n.yub, iter=1:length(mult_x))
        # if length(mult_x) > 10
        # 	break
        # end
    end 
    println(m)
    #solve m again to obtain the optimal fist stage solution and lower bound
    solve(m)
    n.x_opt = getvalue(getindex(m, :x))
    n.y_opt = getvalue(getindex(m, :y))
    n.v_node = getobjectivevalue(m)

    #update multipliers
    n.mult_x = mult_x
    n.mult_y = mult_y
    n.g = g 

    #update upper bound 
    ubsub_problem= []
    for s in scenarios
    	push!(ubsub_problem, generate_ubsub(xbar=n.x_opt, ybar=n.y_opt, Carea=Carea[:,s], prob = prob[s]))
    	a = now()
    	solve(ubsub_problem[s])
    	b = now()
    	ubsub_time = ubsub_time + b - a 
    end

    V_node = getvalue(getindex(m, :firststagecost))
    for s in scenarios
        V_node = V_node + getobjectivevalue(ubsub_problem[s])
    end
    n.V_node = V_node

    n.ubsub_time = ubsub_time
    n.sub_time = sub_time
    n.master_time = master_time
end

#return the most fractional variable's fraction in an array and its index 
function get_frac_var(var_array)
    most_frac_index = 1
    frac = 0.0
    index = 1
    for var in var_array
        if var > frac
            most_frac_index = index
            frac = var 
        end
        index = index + 1
    end
    d = Dict()
    d[:frac] = frac
    d[:index] = most_frac_index
    return d
end























