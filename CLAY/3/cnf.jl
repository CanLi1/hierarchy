using JuMP
using Ipopt
using Mosek
function generate_cnf(;xbar=zeros(length(rectangles)), ybar=zeros(length(rectangles)), Carea=zeros(areas), prob=0)
	# s1 = Model(solver=KnitroSolver())
    # s1 = Model(solver=IpoptSolver())
    s1 = Model(solver=MosekSolver(MSK_IPAR_NUM_THREADS=1, MSK_IPAR_INTPNT_MAX_ITERATIONS=1000000))
   	
   	ϵ=1e-6
   	disjunction = rectangles
   	max_djc = areas
   	@variable(s1, x[i in rectangles])
 	@variable(s1, y[i in rectangles])

 	#second stage var 
 	@variable(s1, xs[i in rectangles])
 	@variable(s1, ys[i in rectangles])
 	@variable(s1, 0<=W[d in disjunction, k in max_djc]<=1)
 	@variable(s1, rsqr[a in areas]>=0)

 	@variable(s1, dot_x[i in rectangles, d in disjunction, k in max_djc])
 	@variable(s1, dot_y[i in rectangles, d in disjunction, k in max_djc])
 	@variable(s1, dot_xs[i in rectangles, d in disjunction, k in max_djc])
 	@variable(s1, dot_ys[i in rectangles, d in disjunction, k in max_djc])
 	@variable(s1, dot_rsqr[a in areas, d in disjunction, k in max_djc]>=0)

 	#x=x̂
 	@constraint(s1, t1[i in rectangles], x[i] == xbar[i])
 	@constraint(s1, t2[i in rectangles], y[i] == ybar[i])


 	# Ax+g(y) ⩽0, g2(y)⩽0
 	@constraint(s1, c1[i in rectangles, j in rectangles, d in disjunction, k in max_djc; i<j], dot_x[i, d, k] - dot_x[j,d,k] == dot_xs[i,d,k]-dot_xs[j,d,k])
 	@constraint(s1, c2[i in rectangles, j in rectangles, d in disjunction, k in max_djc; i<j], dot_y[i, d, k] - dot_y[j,d,k] == dot_ys[i,d,k]-dot_ys[j,d,k])
 	@NLconstraint(s1, c3[d in disjunction, k in max_djc], ((1-ϵ)*W[d,k] + ϵ)*(dot_xs[d,d,k]/((1-ϵ)*W[d,k] + ϵ) - L[d]/2 - xbar[k])^2 - ϵ  * (-L[d]/2-xbar[k])^2*(1-W[d,k]) + ((1-ϵ)*W[d,k] + ϵ) * (dot_ys[d,d,k]/((1-ϵ)*W[d,k] + ϵ) + H[d]/2 - ybar[k])^2 - ϵ *(H[d]/2 - ybar[k])^2*(1-W[d,k]) <= dot_rsqr[k,d,k] )
 	@NLconstraint(s1, c4[d in disjunction, k in max_djc], ((1-ϵ)*W[d,k] + ϵ)*(dot_xs[d,d,k]/((1-ϵ)*W[d,k] + ϵ) - L[d]/2 - xbar[k])^2 - ϵ  * (-L[d]/2-xbar[k])^2*(1-W[d,k]) + ((1-ϵ)*W[d,k] + ϵ) * (dot_ys[d,d,k]/((1-ϵ)*W[d,k] + ϵ) - H[d]/2 - ybar[k])^2 - ϵ *(-H[d]/2 - ybar[k])^2*(1-W[d,k]) <= dot_rsqr[k,d,k] )
 	@NLconstraint(s1, c5[d in disjunction, k in max_djc], ((1-ϵ)*W[d,k] + ϵ)*(dot_xs[d,d,k]/((1-ϵ)*W[d,k] + ϵ) + L[d]/2 - xbar[k])^2 - ϵ  * (+L[d]/2-xbar[k])^2*(1-W[d,k]) + ((1-ϵ)*W[d,k] + ϵ) * (dot_ys[d,d,k]/((1-ϵ)*W[d,k] + ϵ) + H[d]/2 - ybar[k])^2 - ϵ *(H[d]/2 - ybar[k])^2*(1-W[d,k]) <= dot_rsqr[k,d,k] )
 	@NLconstraint(s1, c6[d in disjunction, k in max_djc], ((1-ϵ)*W[d,k] + ϵ)*(dot_xs[d,d,k]/((1-ϵ)*W[d,k] + ϵ) + L[d]/2 - xbar[k])^2 - ϵ  * (+L[d]/2-xbar[k])^2*(1-W[d,k]) + ((1-ϵ)*W[d,k] + ϵ) * (dot_ys[d,d,k]/((1-ϵ)*W[d,k] + ϵ) - H[d]/2 - ybar[k])^2 - ϵ *(-H[d]/2 - ybar[k])^2*(1-W[d,k]) <= dot_rsqr[k,d,k] )

 	#x=∑ẋ, y=∑ẏ, ∑λ=1
 	@constraint(s1, e1[i in rectangles, d in disjunction], sum(dot_x[i,d,k] for k in max_djc)== x[i])
 	@constraint(s1, e2[i in rectangles, d in disjunction], sum(dot_y[i,d,k] for k in max_djc) == y[i])
 	@constraint(s1, e3[i in rectangles, d in disjunction], sum(dot_xs[i,d,k] for k in max_djc) == xs[i])
 	@constraint(s1, e4[i in rectangles, d in disjunction], sum(dot_ys[i,d,k] for k in max_djc) == ys[i])
 	@constraint(s1, e5[a in areas, d in disjunction], sum(dot_rsqr[a,d,k] for k in max_djc) == rsqr[a])
 	@constraint(s1, e6[d in disjunction], sum(W[d,k] for k in max_djc) == 1)

 	#0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
 	@constraint(s1, b1[i in rectangles, d in disjunction, k in max_djc], dot_x[i,d,k]<= 10 * W[d,k])
 	@constraint(s1, b2[i in rectangles, d in disjunction, k in max_djc], dot_x[i,d,k] >= -10 * W[d,k])
 	@constraint(s1, b3[i in rectangles, d in disjunction, k in max_djc], dot_y[i,d,k]<= 40 * W[d,k])
 	@constraint(s1, b4[i in rectangles, d in disjunction, k in max_djc], dot_y[i,d,k] >= -5 * W[d,k])
 	@constraint(s1, b5[i in rectangles, d in disjunction, k in max_djc], dot_xs[i,d,k]<= 10 * W[d,k])
 	@constraint(s1, b6[i in rectangles, d in disjunction, k in max_djc], dot_xs[i,d,k] >= -10 * W[d,k])
 	@constraint(s1, b7[i in rectangles, d in disjunction, k in max_djc], dot_ys[i,d,k]<= 40 * W[d,k])
 	@constraint(s1, b8[i in rectangles, d in disjunction, k in max_djc], dot_ys[i,d,k] >= -5 * W[d,k])
 	@constraint(s1, b9[a in areas, d in disjunction, k in max_djc], dot_rsqr[a,d,k] <= 500*W[d,k])

 	@objective(s1, Min,  prob * sum(rsqr[a] * Carea[a] for a in areas))

 	return s1 
 end









