using JuMP
using Ipopt
function generate_model()
	model = Model(solver=IpoptSolver())
	I = 1:3
	@variable(model, x[i in I]>=1)
	@objective(model, Min, sum(x[i] for i in I ))
	return model
end


function psolve(m::JuMP.Model)
	temp = solve(m)
	d = Dict()
	d[:status] = temp
	d[:x] = getvalue(getindex(m, :x))
	return d
end