addprocs(2)
@everywhere include("input.jl")
@everywhere include("sub.jl")
@everywhere include("master.jl")
@everywhere include("ubsub.jl")
@everywhere include("nlprelax.jl")
@everywhere include("util.jl")
@everywhere include("mosek_sub.jl")
djc_scenarios = []
scenarios = 1:2
for s in scenarios
    push!(djc_scenarios, [[1],[2],[3],[4],[5],[6],[7],[8]])
end

sub_problem=[]
for s in scenarios
    push!(sub_problem, generate_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi, prob=prob[s]))


end
solve(sub_problem[1])
for p in plant
    for i in process
        setvalue(getindex(sub_problem[1], :Qbar)[p,i], -1.0)
        setvalue(getindex(sub_problem[1], :xbar)[p,i], 0)
    end
end

for p in plant
    for i in process
        setvalue(getindex(sub_problem[2], :Qbar)[p,i], 0.0)
        setvalue(getindex(sub_problem[2], :xbar)[p,i], 0)
    end
end
solve(sub_problem[2])
solve(sub_problem[1])

global results = pmap(psolve_sub, sub_problem)

for s in scenarios
    if !results[s][:error]
        sub_problem[s] = results[s][:model]
    end
    # if results[s][:status] != :Optimal
    #     error("NLP solver converges to an infeasible solution")
    # end  
    if results[s][:status] != :Optimal
        sub_problem[s] = generate_mosek_sub(djc=djc_scenarios[s], demand=D[1:2, 1:6, s], price=phi, prob=prob[s])
        temp = solve(sub_problem[s])
        results[s][:status] = temp
        results[s][:QE_dual] = getdual(getindex(sub_problem[s], :t1))
        results[s][:x_dual] = getdual(getindex(sub_problem[s], :t2))
        results[s][:objective] = getobjectivevalue(sub_problem[s])
    end

    
end