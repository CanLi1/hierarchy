include("input.jl")
include("subb.jl")
include("master.jl")
include("ubsub.jl")
include("nlprelax.jl")

#generate subproblem
sub_problem = []
ub_problem = []
sub_stat = []
ub_stat = []
UB = 1e6
LB = -1e6
djc_scenarios = []
mult_QE = []
mult_x = []
g = []
xbar_record = []
QEbar_record = []
time_subproblem = []
for s in scenarios
    push!(djc_scenarios, [[1],[2],[3],[4],[5],[6],[7],[8]])
end

#generate subproblem and master problem 
for s in scenarios
    push!(sub_problem, generate_sub(demand=D[1:2, 1:6, s], price=phi, prob=prob[s]))
end

m = generate_master()
obj_master = []
while UB > LB + 1e-3
    solve(m)
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
    end

    for s in scenarios
    	#update first stage decisions

        a = now()
        for p in plant
            for i in process
                setvalue(getindex(sub_problem[s], :Qbar)[p,i], QEbarr[p,i])
                setvalue(getindex(sub_problem[s], :xbar)[p,i], xbarr[p,i])
            end
        end

        temp = solve(sub_problem[s])

        b = now()

        push!(time_subproblem, b-a)
        if temp != :Optimal
            println(ub_record)
            println(lb_record)
            println(sub_stat)
            println(m)
            println(qbarr)
            println(ybarr)                    
            error("Mosek converge to an infeasible solution")
        end  

        push!(temp_sub_stat, temp)
        QE_dual = getdual(getindex(sub_problem[s], :t1))
        x_dual = getdual(getindex(sub_problem[s], :t2))  	
        for p in plant
        	for i in process
        		temp_mult_QE[s][p,i] = QE_dual[p,i]
        		temp_mult_x[s][p,i] = x_dual[p,i]
        	end
        end
        temp_g[s] = getobjectivevalue(sub_problem[s]) - sum(temp_mult_QE[s][p,i] * QEbarr[p,i] + temp_mult_x[s][p,i] * xbarr[p,i] for p in plant for i in process)
    end

    #update multipliers
    push!(mult_QE, temp_mult_QE)
    push!(mult_x, temp_mult_x)
    push!(g, temp_g)
    push!(sub_stat, temp_sub_stat)

    m = generate_master(mult_QE=mult_QE, mult_x=mult_x, g=g, iter=1:length(mult_QE))
    if length(mult_QE)> 10 #&& obj_master[length(mult_QE)] - obj_master[length(mult_QE)-1] <1e-3
    	break
    end
end 

println(obj_master)
println(time_subproblem)
println(m)
println(mult_QE)
println(mult_x)
println(g)
# println(sub_stat)
# println(mult_QE)
# println(mult_x)
# println(g)
println(xbar_record)
println(QEbar_record)














