include("input.jl")
include("fullspace.jl")

m = generate_fullspace()
solve(m, relaxation=true)
println(getobjectivevalue(m))
println(getvalue(getindex(m, :QE)))
println(getvalue(getindex(m, :x)))
# println(getvalue(getindex(m, :Slack)))
# println(getvalue(getindex(m, :PU)))
# println(getvalue(getindex(m, :F)))
# println(getvalue(getindex(m, :WW)))
# println(getvalue(getindex(m, :y)))
# println(getvalue(getindex(m, :z)))
