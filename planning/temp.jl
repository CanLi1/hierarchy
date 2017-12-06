using JuMP

include("sub.jl")

q = [1.0 2.0 3.0 4.0]
q[1] = 1.9052079022873567
q[2] = 2.2171309919043067
q[3] = 1.5864315019402662
q[4] = 3.200000000000001
y = [1.0 1.0 1.0 1.0]


s = generate_sub(djc=[[1], [2], [3], [4]], demand=6, qbar=q, ybar=y, prob=0.3)
temp = solve(s)
println(temp)
println(getdual(getindex(s, :tq)))