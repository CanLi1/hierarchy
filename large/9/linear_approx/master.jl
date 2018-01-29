using JuMP
# using CPLEX

function generate_master(; mult_x=[], mult_QE =[], g=[], iter=[] )
	# m1 = Model(solver=CplexSolver())
	m1 = Model(solver=GurobiSolver())
	@variable(m1, QE[p in plant, i in process]>=0.0)
	@variable(m1, x[p in plant, i in process], Bin)

	@variable(m1, η[s in scenarios]>=0)

	@constraint(m1, e1[p in plant, i in process], QE[p,i] <= QEU[p,i] * x[p,i])

	#Benders cuts
	@constraint(m1, b1[s in scenarios, it in iter], η[s] >= g[it][s] + sum(mult_x[it][s][p,i] * x[p,i] + mult_QE[it][s][p,i] * QE[p,i] for p in plant for i in process))

	#obj
	@objective(m1, Min, sum(η[s] for s in scenarios))

	return m1
end