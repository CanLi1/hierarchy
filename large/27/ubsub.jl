using Pajarito 
using CPLEX 
# using Ipopt 
# using Mosek
# using BARON
function generate_ubsub(; demand=zeros(2,6), price=zeros(2,6), prob=0.0, Qbar=zeros(2,4), xbar=zeros(2,4))
	#set up model
    s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, log_level=0, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0 ), cont_solver=IpoptSolver(print_level=0)))
    # s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0), cont_solver=MosekSolver()))
    # s1 = Model(solver=KnitroSolver())
    # s1 = Model(solver=BaronSolver(Threads=24))

    @variable(s1, PU[r in supplier, p in plant, j in chemical; (r,j) in RJ]>=0.0)
	@variable(s1, F[p in plant, c in customer, j in chemical]>=0.0)
	@variable(s1, QE[p in plant, i in process]>=0.0)
	@variable(s1, theta[p in plant, i in process, j in chemical, s in scheme]>=0.0)
	@variable(s1, WW[p in plant, i in process, j in chemical, s in scheme]>=0.0)
	@variable(s1, Slack[c in customer, j in chemical]>=0.0)
	@variable(s1, x[p in plant, i in process], Bin)
	@variable(s1, y[r in supplier, p in plant], Bin)
	@variable(s1, z[p in plant, c in customer], Bin)


    # x=x̄
    @constraint(s1, t1[p in plant, i in process], QE[p,i]== Qbar[p,i])
    @constraint(s1, t2[p in plant, i in process], x[p,i]== xbar[p,i])

	@constraint(s1, e1[p in plant, i in process], QE[p,i] <= QEU[p,i] * x[p,i])
	@constraint(s1, e3[p in plant, j in chemical], sum( PU[r,p,j] for r in supplier if (r,j) in RJ) + sum(WW[p,i,j,s] for i in process, s in scheme if ( (i,j) in OJ && (i,s) in PS) ) == sum(F[p,c,j] for c in customer) + sum(WW[p,i,j,s]  for i in process, s in scheme if ((i,j) in IJ && (i,s) in PS)))
	@constraint(s1, e4[p in plant, i in process], sum(theta[p,i,j,s] for j in chemical, s in scheme if ((i,s,j) in JM && (i,s) in PS)) <= H * QE[p,i] * 100.0 )
	@constraint(s1, e5[p in plant, i in process, j in chemical, s in scheme; (i,s) in PS && (i,s,j) in JM], WW[p,i,j,s] == rho * theta[p,i,j,s])
	@constraint(s1, e6[p in plant, i in process, j in chemical, s in scheme; !((i,s) in PS && (i,s,j) in JM)], theta[p,i,j,s] == 0.0)
	@constraint(s1, e7[p in plant, i in process, j in chemical, jj in chemical, s in scheme; (i,s,j) in L && (i,s) in PS && (i,s,jj) in JM], WW[p,i,j,s] == mu[i,s,j] * WW[p,i,jj,s])
	@NLconstraint(s1, e8[p in plant, i in process, j in chemical, jj in chemical, s in scheme; (i,s,j) in Lbar && (i,s) in PS && (i,s,jj) in JM], log(1.0+WW[p,i,j,s]) >= mu[i,s,j] * WW[p,i,jj,s])
	@constraint(s1, e9[r in supplier, p in plant, j in chemical; (r,j) in RJ], PU[r,p,j] <= PUU * y[r,p])
	@constraint(s1, e10[p in plant, c in customer, j in chemical], F[p,c,j] <= FUU * z[p,c])
	@constraint(s1, e11[c in customer, j in chemical], sum(F[p,c,j] for p in plant) + Slack[c,j] == demand[c,j])
	@constraint(s1, e12[p in plant, i in process, j in chemical, s in scheme; !((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)], WW[p,i,j,s] == 0.0)

	@objective(s1, Min, prob * (sum(betaC[i] * QE[p,i] * 100.0 + alphaC[i] * x[p,i] for p in plant, i in process) + sum(delta[i,s] * rho * theta[p,i,j,s] for p in plant, i in process, s in scheme, j in chemical if ((i,s) in PS && (i,s,j) in JM)) + sum((betaS[r,j] + betaRP[r,p])  * PU[r,p,j] for p in plant, j in chemical, r in supplier if (r,j) in RJ) + sum(alphaRP[r,p] * y[r,p] for r in supplier, p in plant) + sum(alphaPC[p,c] * z[p,c] for p in plant, c in customer) + sum(betaPC[p,c] * F[p,c,j] for p in plant, c in customer, j in chemical) + sum(price[c,j] * Slack[c,j] for c in customer, j in chemical) ))
	
	return s1
end