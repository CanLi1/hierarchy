
function generate_nlprelax(; yfbar=zeros(length(integer), length(stages)), nbar=zeros(length(stages)), vbar=log(VL) * ones(length(stages)), Q=zeros(5), prob=0.0)
	# s1 = Model(solver=GurobiSolver(Method=1))
	s1 = Model(solver=CplexSolver(CPX_PARAM_THREADS=1))


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
    @variable(s1, tl_b[i in products])
    @variable(s1, n_v[j in stages])


	#x=x̂
	@constraint(s1, t1[int in integer, j in stages], yf[int, j] == yfbar[int, j])
	@constraint(s1, t2[j in stages], v[j] == vbar[j])
	@constraint(s1, t3[j in stages], n[j] == nbar[j])

	#equations
	@constraint(s1, intp1[i in products, int in intp], tl_b[i]>= exp(int_tl_b[i, int]) + exp(int_tl_b[i, int]) *(tl[i] - b[i] - int_tl_b[i, int])   )
	@constraint(s1, intp2[j in stages, int in intp], n_v[j] >= exp(int_n_v[j, int]) + exp(int_n_v[j, int]) * (n[j] + beta[j] * v[j] - int_n_v[j, int]) )	
	@constraint(s1, e1[j in stages], sum(yf[int,j] for int in integer)== 1 )
	@constraint(s1, e2[j in stages], n[j] == sum(log(int) * yf[int,j] for int in integer))
	@constraint(s1, e3[i in products, j in stages], v[j] >= log(S[i,j]) + b[i])
	@constraint(s1, e4[j in stages], ns[j] == sum(log(int) * ys[int,j] for int in integer))
	@constraint(s1, e5[j in stages], ns[j] <= n[j])
	@constraint(s1, e6[i in products, j in stages], ns[j] + tl[i] >= log(t[i,j]))
	@constraint(s1, e7, sum(Q[i] * tl_b[i]  for i in products) <= H + L)
	@constraint(s1, e8[j in stages], sum(ys[int,j] for int in integer) == 1)

	@constraint(s1, e9, firststagecost>= prob * sum(alpha[j] * n_v[j] for j in stages))
	@constraint(s1, e10, secondstagecost>=prob* ( sum( lambda[j] * sum(k * ys[k, j] for k in integer) for j in stages)+ delta*L) )

	@objective(s1, Min, firststagecost + secondstagecost )	

	return s1 
end