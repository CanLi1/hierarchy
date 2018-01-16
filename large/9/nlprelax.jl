
# function generate_nlprelax(;  demand=zeros(2,6), price=zeros(2,6), xbar=zeros(2,4), Qbar=zeros(2,4), prob=0.0)
function generate_nlprelax(;  demand=zeros(2,6), price=zeros(2,6),  prob=0.0)
	# s1 = Model(solver=MosekSolver(MSK_IPAR_INTPNT_MAX_ITERATIONS=1000000, MSK_DPAR_INTPNT_TOL_REL_GAP=1e-5, MSK_DPAR_INTPNT_NL_TOL_REL_GAP=1e-5))
	s1 = Model(solver=IpoptSolver(print_level=0))

	@NLparameter(s1, xbar[p in plant, i in process]==0.0)
	@NLparameter(s1, Qbar[p in plant, i in process]==0.0)

	s1[:xbar] = xbar
	s1[:Qbar] = Qbar

    @variable(s1, PU[r in supplier, p in plant, j in chemical; (r,j) in RJ]>=0.0)
	@variable(s1, F[p in plant, c in customer, j in chemical]>=0.0)
	@variable(s1, QE[p in plant, i in process]>=0.0)
	@variable(s1, theta[p in plant, i in process, j in chemical, s in scheme; ((i,s) in PS && (i,s,j) in JM)]>=0.0)
	@variable(s1, WW[p in plant, i in process, j in chemical, s in scheme; ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)]>=0.0)
	@variable(s1, Slack[c in customer, j in chemical]>=0.0)
	@variable(s1, 0.0<=x[p in plant, i in process]<=1.0)
	@variable(s1, 0<=y[r in supplier, p in plant]<=1.0)
	@variable(s1, 0<=z[p in plant, c in customer]<=1.0)
	@variable(s1, 0.0<=xx[p in plant, i in process]<=1.0)
	@variable(s1, QEE[p in plant, i in process]>=0.0)
    # x=x̄
	@NLconstraint(s1, t01[p in plant, i in process], QEE[p,i]== Qbar[p,i])
	@NLconstraint(s1, t02[p in plant, i in process], xx[p,i]== xbar[p,i])
	@constraint(s1, t1[p in plant, i in process], QE[p,i]== QEE[p,i])
	@constraint(s1, t2[p in plant, i in process], x[p,i]== xx[p,i])

	#original second stage constraints	
	@constraint(s1, e1[p in plant, i in process], QE[p,i] <= QEU[p,i] * x[p,i])
	@constraint(s1, e3[p in plant, j in chemical], sum( PU[r,p,j] for r in supplier if (r,j) in RJ) + sum(WW[p,i,j,s] for i in process, s in scheme if ( (i,j) in OJ && (i,s) in PS && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar)) ) == sum(F[p,c,j] for c in customer) + sum(WW[p,i,j,s]  for i in process, s in scheme if ((i,j) in IJ && (i,s) in PS && ((i,s,j) in JM || (i,s,j) in L || (i,s,j) in Lbar) )))
	@constraint(s1, e4[p in plant, i in process], sum(theta[p,i,j,s] for j in chemical, s in scheme if ((i,s,j) in JM && (i,s) in PS)) <= H * QE[p,i] * 100.0 )
	@constraint(s1, e5[p in plant, i in process, j in chemical, s in scheme; (i,s) in PS && (i,s,j) in JM ], WW[p,i,j,s] == rho * theta[p,i,j,s])
	@constraint(s1, e7[p in plant, i in process, j in chemical, jj in chemical, s in scheme; (i,s,j) in L && (i,s) in PS && (i,s,jj) in JM], WW[p,i,j,s] == mu[i,s,j] * WW[p,i,jj,s])
	@NLconstraint(s1, e8[p in plant, i in process, j in chemical, jj in chemical, s in scheme; (i,s,j) in Lbar && (i,s) in PS && (i,s,jj) in JM], log(1.0+WW[p,i,j,s]) >= mu[i,s,j] * WW[p,i,jj,s])
	@constraint(s1, e9[r in supplier, p in plant, j in chemical; (r,j) in RJ], PU[r,p,j] <= PUU * y[r,p])
	@constraint(s1, e10[p in plant, c in customer, j in chemical], F[p,c,j] <= FUU * z[p,c])
	@constraint(s1, e11[c in customer, j in chemical], sum(F[p,c,j] for p in plant) + Slack[c,j] == demand[c,j])

	@objective(s1, Min, prob * (sum(betaC[i] * QE[p,i] * 100.0 + alphaC[i] * x[p,i] for p in plant, i in process) + sum(delta[i,s] * rho * theta[p,i,j,s] for p in plant, i in process, s in scheme, j in chemical if ((i,s) in PS && (i,s,j) in JM)) + sum((betaS[r,j] + betaRP[r,p])  * PU[r,p,j] for p in plant, j in chemical, r in supplier if (r,j) in RJ) + sum(alphaRP[r,p] * y[r,p] for r in supplier, p in plant) + sum(alphaPC[p,c] * z[p,c] for p in plant, c in customer) + sum(betaPC[p,c] * F[p,c,j] for p in plant, c in customer, j in chemical) + sum(price[c,j] * Slack[c,j] for c in customer, j in chemical) ) )

	return s1
end