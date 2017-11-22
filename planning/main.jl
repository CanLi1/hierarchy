using JuMP

include("master.jl")
include("ubsub.jl")
include("sub.jl")

m = generate_master()

prob = [0.3 0.4 0.3]
demand = [10 8.5  6]

#generate subproblem
sub_problem = []
ub_problem = []
sub_stat = []
ub_stat = []
obj_sub = 0 
obj_ub = 0
for s = 1:3
    push!(sub_problem, generate_sub(djc = [[1, 2],[3],[4]], demand=demand[s], prob=prob[s], qbar=[2.0797 2.227 1.90697 2.38774], ybar=[1 1 1 1]))
    temp = solve(sub_problem[s])
    push!(sub_stat, temp)
    obj_sub += getobjectivevalue(sub_problem[s])
end
    

# for s = 1:3
#     push!(ub_problem, generate_ubsub(demand=demand[s], prob=prob[s],qbar=[2.0797 2.227 1.90697 2.38774], ybar=[1 1 1 1]))
#     temp = solve(ub_problem[s])
#     push!(ub_stat, temp)
#     obj_ub += getobjectivevalue(ub_problem[s])
# end

println(obj_sub)
println(sub_stat)
# println(sub_problem[1])
# println(obj_ub)
