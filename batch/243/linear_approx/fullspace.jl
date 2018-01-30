using JuMP
using CPLEX 
using Gurobi
include("input.jl")
# include("sub.jl")
# include("ubsub.jl")
# include("util.jl")
# include("nlprelax.jl")

function generate_fullspace()
	# m = Model(solver=PajaritoSolver(rel_gap=0.00001, timeout=10000, mip_solver=CplexSolver(), cont_solver=IpoptSolver()))
	# m = Model(solver=CplexSolver())
	m = Model(solver=GurobiSolver())
	@variable(m, yf[int in integer, j in stages], Bin)
	@variable(m, ys[int in integer, j in stages, w in scenarios], Bin)
	@variable(m, n[j in stages])
	@variable(m, log(VL)<=v[j in stages]<=log(VU))
	@variable(m, ns[j in stages, w in scenarios])
	@variable(m, tl[i in products, w in scenarios])
	@variable(m, b[i in products, w in scenarios])
	@variable(m, L[w in scenarios]>=0)
	@variable(m, firststagecost)
	@variable(m, secondstagecost)
	@variable(m, tl_b[i in products, w in scenarios])
	@variable(m, n_v[j in stages])
	# @variable(m, ns_int[j in stages, w in scenarios])

	#equations
	@constraint(m, e1[j in stages], sum(yf[int,j] for int in integer)== 1 )
	@constraint(m, e2[j in stages], n[j] == sum(log(int) * yf[int,j] for int in integer))
	@constraint(m, e3[i in products, j in stages, w in scenarios], v[j] >= log(S[i,j]) + b[i,w])
	@constraint(m, e4[j in stages, w in scenarios], ns[j,w] == sum(log(int) * ys[int,j,w] for int in integer))
	@constraint(m, e5[j in stages, w in scenarios], ns[j,w] <= n[j])
	@constraint(m, e6[i in products, j in stages, w in scenarios], ns[j,w] + tl[i,w] >= log(t[i,j]))

	#constraints for intp 
	@constraint(m, intp1[i in products, w in scenarios, int in intp], tl_b[i, w]>= exp(int_tl_b[i, int]) + exp(int_tl_b[i, int]) *(tl[i,w] - b[i,w] - int_tl_b[i, int])   )
	@constraint(m, intp2[j in stages, int in intp], n_v[j] >= exp(int_n_v[j, int]) + exp(int_n_v[j, int]) * (n[j] + beta[j] * v[j] - int_n_v[j, int]) )
	# @constraint(m, intp3[j in stages, int in integer, w in scenarios], ns_int[j,w] >=  exp(int_ns[int, j, w]) + exp(int_ns[int, j, w]) * (ns[j,w] - int_ns[int, j, w])  )
	# @NLconstraint(m, e7[w in scenarios], sum(Q[i,w] * exp(tl[i,w]-b[i,w]) for i in products) <= H + L[w])
	@constraint(m, e7[w in scenarios], sum(Q[i,w] * tl_b[i,w]  for i in products) <= H + L[w])
	@constraint(m, e8[j in stages, w in scenarios], sum(ys[int,j,w] for int in integer) == 1)

	@constraint(m, e9, firststagecost>= sum(alpha[j] * n_v[j] for j in stages))
	@constraint(m, e10, secondstagecost>=sum(prob[w] * ( sum( lambda[j] * sum(k * ys[k, j,w] for k in integer) for j in stages)+ delta*L[w]) for w in scenarios ))

	@objective(m, Min, firststagecost + secondstagecost )
	return m 
end

m = generate_fullspace()
solve(m)
yf = getvalue(getindex(m, :yf))
ys = getvalue(getindex(m, :ys))
v = getvalue(getindex(m, :v))
L = getvalue(getindex(m, :L))
# println(yf)
# println(v)
# println(L)
# println(ys)
println(getobjectivevalue(m))
# yf = getvalue(getindex(m, :yf))
# n = zeros(6)
# for j in stages
# 	n[j] = sum(log(k) * yf[k,j] for k in integer)
# end
# v = getvalue(getindex(m, :v))
# # v = ones(6) * log(VL)
# # n = ones(6) * log(1)
# # yf = zeros(4,6)
# # for j in stages
# # 	yf[1,j] = 1
# # end

# relax = 0.0
# ub = 0.0
# relax_sub = []
# ub_sub_problem = []
# f_cost = 0.0 
# s_cost = 0.0 
# L_record = []
# ns_record = []
# firts_cost_scenario = []
# second_cost_scenario = []
# for s in scenarios
# 	push!(relax_sub, generate_sub(Q=Q[:,s], prob =prob[s]))
# 	push!(ub_sub_problem, generate_ubsub(Q=Q[:,s], prob=prob[s], nbar=n, vbar=v, yfbar=yf))
# 	for j in stages
# 		setvalue(getindex(relax_sub[s], :nbar)[j], n[j])
# 		setvalue(getindex(relax_sub[s], :vbar)[j], v[j])		
# 		for int in integer
# 			setvalue(getindex(relax_sub[s], :yfbar)[int, j], yf[int, j])
# 		end
# 	end
# 	solve(relax_sub[s])
# 	solve(ub_sub_problem[s])
# 	relax = relax + getobjectivevalue(relax_sub[s])
# 	ub = ub + getobjectivevalue(ub_sub_problem[s])
# 	f_cost = f_cost + getvalue(getindex(ub_sub_problem[s], :firststagecost))
# 	s_cost = s_cost + getvalue(getindex(ub_sub_problem[s], :secondstagecost))
# 	push!(firts_cost_scenario, getvalue(getindex(ub_sub_problem[s], :firststagecost)))
# 	push!(second_cost_scenario, getvalue(getindex(ub_sub_problem[s], :secondstagecost)))
# 	push!(L_record, getvalue(getindex(ub_sub_problem[s], :L)))
# 	push!(ns_record, getvalue(getindex(ub_sub_problem[s], :ns)))

# end

# println(relax)
# println(ub)
# # println(getvalue(getindex(m, :firststagecost)))
# # println(getvalue(getindex(m, :secondstagecost)))
# println(f_cost)
# println(s_cost)
# println(firts_cost_scenario)
# println(second_cost_scenario)
# # println(getvalue(getindex(m, :L)))
# # println(getvalue(getindex(m, :ns)))
# println(L_record)
# println(ns_record)



















