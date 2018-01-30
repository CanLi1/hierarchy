function generate_nlprelax(; Q=zeros(5), prob=0.0)
	s1 = Model(solver=IpoptSolver(print_level=0))
	# s1 = Model(solver=MosekSolver())
	# s1 = Model(solver=BaronSolver())

	@NLparameter(s1, yfbar[int in  integer, j in stages]==0.0)
	@NLparameter(s1, nbar[j in stages]==0.0)
	@NLparameter(s1, vbar[j in stages]==0.0)

	s1[:yfbar] = yfbar
	s1[:nbar] = nbar
	s1[:vbar] = vbar

	@variable(s1, 0<=yf[int in integer, j in stages]<=1)
	@variable(s1, 0<=ys[int in integer, j in stages]<=1)
	@variable(s1, n[j in stages])
	@variable(s1, log(VL)<=v[j in stages]<=log(VU))
	@variable(s1, ns[j in stages])
	@variable(s1, tl[i in products])
	@variable(s1, b[i in products])
	@variable(s1, L>=0)
	@variable(s1, firststagecost)
	@variable(s1, secondstagecost)

	@variable(s1, 0<=yff[int in integer, j in stages]<=1)
	@variable(s1, nn[j in stages])
	@variable(s1, log(VL)<=vv[j in stages]<=log(VU))

	#x=x̂
	@NLconstraint(s1, t01[int in integer, j in stages], yff[int, j]==yfbar[int, j])
	@NLconstraint(s1, t02[j in stages], nn[j]==nbar[j])
	@NLconstraint(s1, t03[j in stages], vv[j]==vbar[j])
	@constraint(s1, t1[int in integer, j in stages], yf[int, j] == yff[int, j])
	@constraint(s1, t2[j in stages], v[j] == vv[j])
	@constraint(s1, t3[j in stages], n[j] == nn[j])

	#equations
	@constraint(s1, e1[j in stages], sum(yf[int,j] for int in integer)== 1 )
	@constraint(s1, e2[j in stages], n[j] == sum(log(int) * yf[int,j] for int in integer))
	@constraint(s1, e3[i in products, j in stages], v[j] >= log(S[i,j]) + b[i])
	@constraint(s1, e4[j in stages], ns[j] == sum(log(int) * ys[int,j] for int in integer))
	@constraint(s1, e5[j in stages], ns[j] <= n[j])
	@constraint(s1, e6[i in products, j in stages], ns[j] + tl[i] >= log(t[i,j]))
	@NLconstraint(s1, e7, sum(Q[i] * exp(tl[i]-b[i]) for i in products) <= H + L)
	@constraint(s1, e8[j in stages], sum(ys[int,j] for int in integer) == 1)

	@NLconstraint(s1, e9, firststagecost>= prob * sum(alpha[j] * exp(n[j] + beta[j] * v[j]) for j in stages))
	@NLconstraint(s1, e10, secondstagecost>=prob *( sum( lambda[j] * exp(ns[j]) for j in stages)+ delta*L) )

	@objective(s1, Min, firststagecost + secondstagecost )	

	return s1 
end