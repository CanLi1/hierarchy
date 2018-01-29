include("util.jl")
include("input1.jl")
include("sub.jl")
s= generate_sub(Q=Q[:,1], prob =prob[1])
for j in stages
	JuMP.setRHS(getindex(s, :t1)[1, j], 1)
end
temp = solve(s)
println(temp)