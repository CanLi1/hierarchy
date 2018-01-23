# using Pajarito 
# using CPLEX 
# using Ipopt 
# using Mosek
using BARON
function generate_ubsub(; Q=zeros(length(products)), prob=0.0, nbar=zeros(length(stages)), vbar=zeros(length(stages)), yfbar=zeros(length(integer), length(stages)))
	# s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0), cont_solver=IpoptSolver(print_level=0)))
	s1 = Model(solver=BaronSolver())

	@variable(s1, yf[int in integer, j in stages], Bin)
	@variable(s1, ys[int in integer, j in stages], Bin)
	@variable(s1, n[j in stages])
	@variable(s1, log(VL)<=v[j in stages]<=log(VU))
	@variable(s1, ns[j in stages])
	@variable(s1, tl[i in products])
	@variable(s1, b[i in products])
	@variable(s1, L>=0)
	@variable(s1, firststagecost)
	@variable(s1, secondstagecost)



	#x=x̂
	@constraint(s1, t1[int in integer, j in stages], yf[int, j] == yfbar[int, j])
	@constraint(s1, t2[j in stages], v[j] == vbar[j])
	@constraint(s1, t3[j in stages], n[j] == nbar[j])

	#equations
	@constraint(s1, e3[i in products, j in stages], v[j] >= log(S[i,j]) + b[i])
	@constraint(s1, e4[j in stages], ns[j] == sum(log(int) * ys[int,j] for int in integer))
	@constraint(s1, e5[j in stages], ns[j] <= n[j])
	@constraint(s1, e6[i in products, j in stages], ns[j] + tl[i] >= log(t[i,j]))
	@NLconstraint(s1, e7, sum(Q[i] * exp(tl[i]-b[i]) for i in products) <= H + L)
	@constraint(s1, e8[j in stages], sum(ys[int,j] for int in integer) == 1)

	@NLconstraint(s1, e9, firststagecost>= prob * sum(alpha[j] * exp(n[j] + beta[j] * v[j]) for j in stages))
	@NLconstraint(s1, e10, secondstagecost>=prob * ( sum( lambda[j] * exp(ns[j]) for j in stages)+ delta*L) )

	@NLobjective(s1, Min, firststagecost + secondstagecost )	

	return s1 
end