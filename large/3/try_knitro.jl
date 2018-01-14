include("input.jl")
include("sub.jl")
sub_problem = []
for s in scenarios
    push!(sub_problem, generate_sub(demand=D[1:2, 1:6, s], price=phi, prob=prob[s]))
    solve(sub_problem[s])
end