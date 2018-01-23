using JuMP
using KNITRO
function generate_fullspace()
	# m = Model(solver=PajaritoSolver(rel_gap=0.001, timeout=10000, mip_solver=CplexSolver(), cont_solver=MosekSolver(MSK_IPAR_NUM_THREADS=1)))
	#m = Model(solver=BaronSolver(maxtime=5e4))
	m=Model(solver=KnitroSolver())
	@variable(m, PU[r in supplier, p in plant, j in chemical, w in scenarios; (r,j) in RJ]>=0.0)
	@variable(m, F[p in plant, c in customer, j in chemical, w in scenarios]>=0.0)
	@variable(m, QE[p in plant, i in process]>=0.0)
	@variable(m, theta[p in plant, i in process, j in chemical, s in scheme, w in scenarios; ((i,s) in PS && (i,s,j) in JM)]>=0.0)
	@variable(m, WW[p in plant, i in process, j in chemical, s in scheme, w in scenarios; ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)]>=0.0)
	@variable(m, Slack[c in customer, j in chemical, w in scenarios]>=0.0)
	@variable(m, x[p in plant, i in process], Bin)
	@variable(m, y[r in supplier, p in plant, w in scenarios], Bin)
	@variable(m, z[p in plant, c in customer, w in scenarios], Bin)

	@constraint(m, e1[p in plant, i in process], QE[p,i] <= QEU[p,i] * x[p,i])
	@constraint(m, e3[p in plant, j in chemical, w in scenarios], sum( PU[r,p,j,w] for r in supplier if (r,j) in RJ) + sum(WW[p,i,j,s,w] for i in process, s in scheme if ( (i,j) in OJ && (i,s) in PS && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)) ) == sum(F[p,c,j,w] for c in customer) + sum(WW[p,i,j,s,w]  for i in process, s in scheme if ((i,j) in IJ && (i,s) in PS && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar) )))
	@constraint(m, e4[p in plant, i in process, w in scenarios], sum(theta[p,i,j,s,w] for j in chemical, s in scheme if ((i,s,j) in JM && (i,s) in PS)) <= H * QE[p,i] * 100.0 )
	@constraint(m, e5[p in plant, i in process, j in chemical, s in scheme, w in scenarios; (i,s) in PS && (i,s,j) in JM ], WW[p,i,j,s,w] == rho * theta[p,i,j,s,w])
	# @constraint(m, e6[p in plant, i in process, j in chemical, w in scenarios, s in scheme; !((i,s) in PS && (i,s,j) in JM)], theta[p,i,j,s,w] == 0.0)
	@constraint(m, e7[p in plant, i in process, j in chemical, jj in chemical, w in scenarios, s in scheme; (i,s,j) in L && (i,s) in PS && (i,s,jj) in JM], WW[p,i,j,s,w] == mu[i,s,j] * WW[p,i,jj,s,w])
	@NLconstraint(m, e8[p in plant, i in process, j in chemical, jj in chemical, w in scenarios, s in scheme; (i,s,j) in Lbar && (i,s) in PS && (i,s,jj) in JM], log(1.0+WW[p,i,j,s,w]) >= mu[i,s,j] * WW[p,i,jj,s,w])
	@constraint(m, e9[r in supplier, p in plant, w in scenarios, j in chemical; (r,j) in RJ], PU[r,p,j,w] <= PUU * y[r,p,w])
	@constraint(m, e10[p in plant, c in customer, j in chemical, w in scenarios], F[p,c,j,w] <= FUU * z[p,c,w])
	@constraint(m, e11[c in customer, j in chemical, w in scenarios], sum(F[p,c,j,w] for p in plant) + Slack[c,j,w] == D[c,j,w])
	# @constraint(m, e12[p in plant, i in process, j in chemical, s in scheme, w in scenarios; !((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)], WW[p,i,j,s,w] == 0.0)

	@objective(m, Min, sum(prob[w] * (sum(betaC[i] * QE[p,i] * 100.0 + alphaC[i] * x[p,i] for p in plant, i in process) + sum(delta[i,s] * rho * theta[p,i,j,s,w] for p in plant, i in process, s in scheme, j in chemical if ((i,s) in PS && (i,s,j) in JM)) + sum((betaS[r,j] + betaRP[r,p])  * PU[r,p,j,w] for p in plant, j in chemical, r in supplier if (r,j) in RJ) + sum(alphaRP[r,p] * y[r,p,w] for r in supplier, p in plant) + sum(alphaPC[p,c] * z[p,c,w] for p in plant, c in customer) + sum(betaPC[p,c] * F[p,c,j,w] for p in plant, c in customer, j in chemical) + sum(phi[c,j] * Slack[c,j,w] for c in customer, j in chemical) ) for w in scenarios))
	return m
end		

include("input.jl")

m = generate_fullspace()
solve(m)
println(getobjectivevalue(m))
println(getvalue(getindex(m, :QE)))
println(getvalue(getindex(m, :x)))
