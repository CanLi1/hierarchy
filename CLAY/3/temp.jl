 using JuMP
 using Pajarito
 using CPLEX 
 using Ipopt
include("input.jl")
include("dnf.jl")
include("ubsub.jl")
s = generate_ubsub()
solve(s)
