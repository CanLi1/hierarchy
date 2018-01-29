function generate_master(; mult_n=[], mult_yf =[], mult_v=[], g=[], iter=[] )
	# m1 = Model(solver=CplexSolver())
	m1 = Model(solver=GurobiSolver())
	@variable(m1, yf[int in integer, j in stages], Bin)
	@variable(m1, n[j in stages])
	@variable(m1, log(VL)<=v[j in stages]<=log(VU))

	@variable(m1, η[s in scenarios]>=0)

	@constraint(m1, e1[j in stages], sum(yf[int,j] for int in integer)== 1 )
	@constraint(m1, e2[j in stages], n[j] == sum(log(int) * yf[int,j] for int in integer))
	

	#Benders cuts
	@constraint(m1, b1[s in scenarios, it in iter], η[s] >= g[it][s] + sum(sum(mult_yf[it][s][int, j] * yf[int,j] for int in integer) + mult_v[it][s][j] * v[j] + mult_n[it][s][j]*n[j] for j in stages))

	#obj
	@objective(m1, Min, sum(η[s] for s in scenarios))

	return m1
end