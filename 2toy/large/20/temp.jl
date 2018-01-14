include("input.jl")
include("sub.jl")
p1=3.8211145495655625
p2=1.9838699412249086
x1=1.0
x2=1.000000000000000
obj = 0.0
for s in scenarios
	m = generate_sub(p1bar=p1, p2bar=p2, x1bar=x1, x2bar=x2, demand=d[s], prob=prob[s])
	solve(m)
	obj = obj + getobjectivevalue(m)
end
println(obj)

