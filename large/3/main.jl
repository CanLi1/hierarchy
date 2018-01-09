

addprocs(2)
@everywhere include("input.jl")
@everywhere include("sub.jl")
@everywhere include("master.jl")
@everywhere include("ubsub.jl")
@everywhere include("nlprelax.jl")
@everywhere include("util.jl")
#generate subproblem
sub_problem = []
ub_problem = []
sub_stat = []
ub_stat = []
UB = 1e6
LB = -1e6
relax_UB = 1e6
relax_LB = -1e6
djc_scenarios = []
mult_QE = []
mult_x = []
g = []
xbar_record = []
QEbar_record = []
sub_obj_record = []

for s in scenarios
    push!(djc_scenarios, [[1],[2],[3],[4],[5],[6],[7],[8]])
end

#generate subproblem and master problem 
for s in scenarios
    push!(sub_problem, generate_nlprelax(demand=D[1:2, 1:6, s], price=phi, prob=prob[s]))
end

m = generate_master()
obj_master = []
a = now()
while relax_UB > relax_LB + 1e-1
    solve(m)
    relax_LB = getobjectivevalue(m)
    push!(obj_master, getobjectivevalue(m))
    QEbarr = getvalue(getindex(m, :QE))
    xbarr = getvalue(getindex(m, :x))
    push!(xbar_record, xbarr)
    push!(QEbar_record, QEbarr)
    temp_mult_QE = []
    temp_mult_x = []
    temp_g = []
    temp_sub_stat = []

    #initialize 
    for s in scenarios
    	push!(temp_mult_QE, zeros(length(plant), length(process)))
    	push!(temp_mult_x, zeros(length(plant), length(process)))
    	push!(temp_g, 0.0)

    	#update first stage decisions
        for p in plant
            for i in process
                setvalue(getindex(sub_problem[s], :Qbar)[p,i], QEbarr[p,i])
                setvalue(getindex(sub_problem[s], :xbar)[p,i], xbarr[p,i])
            end
        end
    end

    results = pmap(psolve, sub_problem)

    for s in scenarios
        sub_problem[s] = results[s][:model]
        if results[s][:status] != :Optimal
            error("NLP solver converges to an infeasible solution")
        end  

        push!(temp_sub_stat, results[s][:status])	
        for p in plant
        	for i in process
        		temp_mult_QE[s][p,i] = results[s][:QE_dual][p,i]
        		temp_mult_x[s][p,i] = results[s][:x_dual][p,i]
        	end
        end
        temp_g[s] = results[s][:objective]- sum(temp_mult_QE[s][p,i] * QEbarr[p,i] + temp_mult_x[s][p,i] * xbarr[p,i] for p in plant for i in process)
    end

    #update upper bound 
    if relax_UB > sum(results[s][:objective] for s in scenarios)
        relax_UB = sum(results[s][:objective] for s in scenarios)
    end

    push!(sub_obj_record, sum(results[s][:objective] for s in scenarios))

    #update multipliers
    push!(mult_QE, temp_mult_QE)
    push!(mult_x, temp_mult_x)
    push!(g, temp_g)
    push!(sub_stat, temp_sub_stat)

    m = generate_master(mult_QE=mult_QE, mult_x=mult_x, g=g, iter=1:length(mult_QE))

    # if length(mult_QE)> 50 && obj_master[length(mult_QE)] - obj_master[length(mult_QE)-1] <1e-3
    # 	break
    # end
end

temp_length = length(mult_QE)



QEbar_for_ub = zeros(length(plant), length(process))
xbar_for_ub = zeros(length(plant), length(process))

while UB > LB + 1e-1
    relax_UB = 1e6
    #generate subproblem and master problem 
    sub_problem=[]
    for s in scenarios
        push!(sub_problem, generate_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi, prob=prob[s]))
    end
    m = generate_master(mult_QE=mult_QE, mult_x=mult_x, g=g, iter=1:length(mult_QE))

    while relax_UB > relax_LB + 1e-1
        solve(m)
        push!(obj_master, getobjectivevalue(m))
        relax_LB = getobjectivevalue(m)
        QEbarr = getvalue(getindex(m, :QE))
        xbarr = getvalue(getindex(m, :x))
        push!(xbar_record, xbarr)
        push!(QEbar_record, QEbarr)
        temp_mult_QE = []
        temp_mult_x = []
        temp_g = []
        temp_sub_stat = []

        #initialize 
        temp_obj = 0.0
        for s in scenarios
            push!(temp_mult_QE, zeros(length(plant), length(process)))
            push!(temp_mult_x, zeros(length(plant), length(process)))
            push!(temp_g, 0.0)
         
            for p in plant
                for i in process
                    setvalue(getindex(sub_problem[s], :Qbar)[p,i], QEbarr[p,i])
                    setvalue(getindex(sub_problem[s], :xbar)[p,i], xbarr[p,i])
                end
            end

        end

        global results = pmap(psolve_sub, sub_problem)

        for s in scenarios
            sub_problem[s] = results[s][:model]
            if results[s][:status] != :Optimal
                error("NLP solver converges to an infeasible solution")
            end  

            push!(temp_sub_stat, results[s][:status])   
            for p in plant
                for i in process
                    temp_mult_QE[s][p,i] = results[s][:QE_dual][p,i]
                    temp_mult_x[s][p,i] = results[s][:x_dual][p,i]
                end
            end
            temp_g[s] = results[s][:objective]- sum(temp_mult_QE[s][p,i] * QEbarr[p,i] + temp_mult_x[s][p,i] * xbarr[p,i] for p in plant for i in process)
        end

        #update upper bound 
        if relax_UB > sum(results[s][:objective] for s in scenarios)
            relax_UB = sum(results[s][:objective] for s in scenarios)
        end
        push!(sub_obj_record, sum(results[s][:objective] for s in scenarios))
        
        #update multipliers
        push!(mult_QE, temp_mult_QE)
        push!(mult_x, temp_mult_x)
        push!(g, temp_g)
        push!(sub_stat, temp_sub_stat)

        m = generate_master(mult_QE=mult_QE, mult_x=mult_x, g=g, iter=1:length(mult_QE))
        if relax_UB > relax_LB + 1e-1
            QEbar_for_ub = QEbarr
            xbar_for_ub = xbarr
        end
        if length(mult_QE) > temp_length + 200
            break
        end
    end 

    #update lower bound
    LB = relax_LB

    #regenerate the sub_problem and solve again
    sub_problem = []
    for s in scenarios
        push!(sub_problem, generate_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi, prob=prob[s]))
        for p in plant
            for i in process
                setvalue(getindex(sub_problem[s], :Qbar)[p,i], QEbar_for_ub[p,i])
                setvalue(getindex(sub_problem[s], :xbar)[p,i], xbar_for_ub[p,i])
            end
        end
        solve(sub_problem[s])
    end    
    #now check if the solution is in convex hull 
    all_scenarios_in_convex_hull = true

    for s in scenarios
        λ = getvalue(getindex(sub_problem[s], :dot_λ))
        dot_y = getvalue(getindex(sub_problem[s], :dot_y))
        cur_djc = djc_scenarios[s]
        index_d = 1
        y_over_λ = []
        is_in_convex_hull = false
        for d in cur_djc
            y_over_λ_d = []
            is_djc_d_in_convex_hull = true
            for k in 1:2^length(d)
                y_over_λ_d_k = -1 * ones(8)
                if λ[index_d, k] > 1e-5
                    for i in J1
                        y_over_λ_d_k[i] = dot_y[i, index_d, k] / λ[index_d, k]
                    end
                end
                if  !check_integer(y_over_λ_d_k)
                    is_djc_d_in_convex_hull = false
                end
                push!(y_over_λ_d, y_over_λ_d_k)
            end
            if is_djc_d_in_convex_hull 
                is_in_convex_hull = true
            end
            push!(y_over_λ, y_over_λ_d)
            index_d = index_d + 1
        end

        #if not in convex hull, apply one basic step
        if !is_in_convex_hull
            all_scenarios_in_convex_hull = false
            #try to merge the first and second disjunction
            temp_djc = []
            for k in 1:(length(cur_djc)-1)
                push!(temp_djc, [0])
            end
            temp_djc[1] = [cur_djc[1]; cur_djc[2]]
            for k in 1:(length(cur_djc)-2)
                temp_djc[k+1] = cur_djc[k+2]
            end 
            djc_scenarios[s] = temp_djc
        end

        println(y_over_λ)
    end
    #obtain an upper bound
    #generate upper bound subproblem
    ub_sub_problem = []
    for s in scenarios
        push!(ub_sub_problem, generate_ubsub(demand=D[1:2, 1:6, s], price=phi, prob=prob[s], Qbar=QEbar_for_ub, xbar=xbar_for_ub))
    end

    #solve upper bound subproblem in parallel
    results = pmap(ub_psolve, ub_sub_problem)
    temp_UB = sum(results[s][:objective] for s in scenarios)
    if UB > temp_UB
        UB = temp_UB
    end
    break

end
 
println(temp_length)
println(UB)
println(LB)
println(obj_master)
println(sub_obj_record)
println(djc_scenarios)

b=now()
println(b-a)
# println(sub_stat)
# println(mult_QE)
# println(mult_x)
# println(g)
# println(xbar_record)
# println(QEbar_record)














