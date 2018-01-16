

addprocs(23)
@everywhere include("input.jl")
@everywhere include("sub.jl")
@everywhere include("master.jl")
@everywhere include("ubsub.jl")
@everywhere include("nlprelax.jl")
@everywhere include("util.jl")
@everywhere include("mosek_sub.jl")
@everywhere include("cnf_result.jl")
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
xbar_record = []
QEbar_record = []
sub_obj_record = []
a = now()
b = now()
sub_time = a - b
master_time = a - b
ubsub_time = a - b 
resolve_sub_time = a - b 

for s in scenarios
    push!(djc_scenarios, [)
end



while UB > LB * 1.001
    relax_UB = 1e6
    #generate subproblem and master problem 
    sub_problem=[]
    for s in scenarios
        push!(sub_problem, generate_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi[1:2, 1:6, s], prob=prob[s]))
    end
    m = generate_master(mult_QE=mult_QE, mult_x=mult_x, g=g, iter=1:length(mult_QE))

    while relax_UB > relax_LB + 1e-1
        c= now()
        solve(m)
        d = now()
        master_time = master_time + d - c 
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
        c = now()
        global results = pmap(psolve_sub, sub_problem)
        d = now()
        sub_time = sub_time + d - c 

        for s in scenarios
            if !results[s][:error] && results[s][:status] == :Optimal
                sub_problem[s] = results[s][:model]
            end
            # if results[s][:status] != :Optimal
            #     error("NLP solver converges to an infeasible solution")
            # end  
            if results[s][:status] != :Optimal
                sub_problem[s] = generate_mosek_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi[1:2, 1:6, s], prob=prob[s])
                #update first stage decisions
                for p in plant
                    for i in process
                        setvalue(getindex(sub_problem[s], :Qbar)[p,i], QEbarr[p,i])
                        setvalue(getindex(sub_problem[s], :xbar)[p,i], xbarr[p,i])
                    end
                end
                c = now()
                temp = solve(sub_problem[s])
                d = now()
                sub_time = sub_time + d - c
                results[s][:status] = temp
                results[s][:QE_dual] = getdual(getindex(sub_problem[s], :t1))
                results[s][:x_dual] = getdual(getindex(sub_problem[s], :t2))
                results[s][:objective] = getobjectivevalue(sub_problem[s])

                #revert back to knitro
                sub_problem[s] = generate_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi[1:2, 1:6, s], prob=prob[s])
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
        push!(sub_problem, generate_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi[1:2, 1:6, s], prob=prob[s]))
        for p in plant
            for i in process
                setvalue(getindex(sub_problem[s], :Qbar)[p,i], QEbar_for_ub[p,i])
                setvalue(getindex(sub_problem[s], :xbar)[p,i], xbar_for_ub[p,i])
            end
        end
        c = now()
        temp = solve(sub_problem[s])
        d = now()
        resolve_sub_time = resolve_sub_time + d - c 
        if temp != :Optimal
            sub_problem[s] = generate_mosek_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi[1:2, 1:6, s], prob=prob[s])
            for p in plant
                for i in process
                    setvalue(getindex(sub_problem[s], :Qbar)[p,i], QEbar_for_ub[p,i])
                    setvalue(getindex(sub_problem[s], :xbar)[p,i], xbar_for_ub[p,i])
                end
            end
            c = now()
            solve(sub_problem[s])
            d = now()
            resolve_sub_time = resolve_sub_time + d - c 
        end

    end    
    
    #obtain an upper bound
    #generate upper bound subproblem
    ub_sub_problem = []
    for s in scenarios
        push!(ub_sub_problem, generate_ubsub(demand=D[1:2, 1:6, s], price=phi[1:2, 1:6, s], prob=prob[s], Qbar=QEbar_for_ub, xbar=xbar_for_ub))
    end

    #solve upper bound subproblem in parallel
    c = now()
    temp_UB = 0.0
    for s in scenarios
        solve(ub_sub_problem[s])
        temp_UB = temp_UB + getobjectivevalue(ub_sub_problem[s])
    end
    d = now()
    ubsub_time = ubsub_time +  d - c 

    if UB > temp_UB
        UB = temp_UB
    end
    


break 
   

end
b=now()
println(b-a)
println(temp_length)
println(UB)
println(LB)
println(obj_master)
println(sub_obj_record)
println(djc_scenarios)
println("Optimal first stage solution")
print("xbar=")
println(xbar_for_ub)
print("Qbar=")
println(QEbar_for_ub)
println("ubsub_time")
println(ubsub_time)
println("sub_time")
println(sub_time)
println("master_time")
println(master_time)
println("resolve_sub_time")
println(resolve_sub_time)

print("mult_QE=")
println(mult_QE)
print("mult_x=")
println(mult_x)
println("g")
println(g)


# println(sub_stat)
# println(mult_QE)
# println(mult_x)
# println(g)
# println(xbar_record)
# println(QEbar_record)














