function generate_knitro_dnf(;xbar=zeros(length(rectangles)), ybar=zeros(length(rectangles)), Carea=zeros(areas), prob=0, xub = zeros(length(rectangles)), yub = zeros(length(rectangles)), xlb = zeros(length(rectangles)), ylb = zeros(length(rectangles)))
	# s1 = Model(solver=KnitroSolver())
    s1 = Model(solver=IpoptSolver())
    # s1 = Model(solver=MosekSolver(MSK_IPAR_NUM_THREADS=1, MSK_IPAR_INTPNT_MAX_ITERATIONS=1000000))
   	
   	ϵ=1e-6
   	max_djc = 1:8
   	@variable(s1, x[i in rectangles])
 	@variable(s1, y[i in rectangles])

 	#second stage var 
 	@variable(s1, xs[i in rectangles])
 	@variable(s1, ys[i in rectangles])
 	@variable(s1, 0<=W[k in max_djc]<=1)
 	@variable(s1, rsqr[a in areas]>=0)

  	@variable(s1, dot_x[i in rectangles,  k in max_djc])
 	@variable(s1, dot_y[i in rectangles,  k in max_djc])
 	@variable(s1, dot_xs[i in rectangles,  k in max_djc])
 	@variable(s1, dot_ys[i in rectangles,  k in max_djc])
 	@variable(s1, dot_rsqr[a in areas,  k in max_djc]>=0)

 	#x=x̂
 	@constraint(s1, t1[i in rectangles], x[i] == xbar[i])
 	@constraint(s1, t2[i in rectangles], y[i] == ybar[i])

  	# Ax+g(y) ⩽0, g2(y)⩽0
 	@constraint(s1, c1[i in rectangles, j in rectangles,  k in max_djc; i<j], dot_x[i,k] - dot_x[j,k] == dot_xs[i,k]-dot_xs[j,k])
 	@constraint(s1, c2[i in rectangles, j in rectangles,  k in max_djc; i<j], dot_y[i,k] - dot_y[j,k] == dot_ys[i,k]-dot_ys[j,k])
 	for i1 in areas
 		for i2 in areas
 			for i3 in areas
 				cur_k = (i1-1)*4 + (i2-1) *2 + i3
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - L[1]/2 - xbar[i1])^2 - ϵ  * (-L[1]/2-xbar[i1])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + H[1]/2 - ybar[i1])^2 - ϵ *(H[1]/2 - ybar[i1])^2*(1-W[cur_k]) <= dot_rsqr[i1,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - L[1]/2 - xbar[i1])^2 - ϵ  * (-L[1]/2-xbar[i1])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - H[1]/2 - ybar[i1])^2 - ϵ *(-H[1]/2 - ybar[i1])^2*(1-W[cur_k]) <= dot_rsqr[i1,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + L[1]/2 - xbar[i1])^2 - ϵ  * (+L[1]/2-xbar[i1])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + H[1]/2 - ybar[i1])^2 - ϵ *(H[1]/2 - ybar[i1])^2*(1-W[cur_k]) <= dot_rsqr[i1,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + L[1]/2 - xbar[i1])^2 - ϵ  * (+L[1]/2-xbar[i1])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[1,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - H[1]/2 - ybar[i1])^2 - ϵ *(-H[1]/2 - ybar[i1])^2*(1-W[cur_k]) <= dot_rsqr[i1,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - L[2]/2 - xbar[i2])^2 - ϵ  * (-L[2]/2-xbar[i2])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + H[2]/2 - ybar[i2])^2 - ϵ *(H[2]/2 - ybar[i2])^2*(1-W[cur_k]) <= dot_rsqr[i2,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - L[2]/2 - xbar[i2])^2 - ϵ  * (-L[2]/2-xbar[i2])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - H[2]/2 - ybar[i2])^2 - ϵ *(-H[2]/2 - ybar[i2])^2*(1-W[cur_k]) <= dot_rsqr[i2,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + L[2]/2 - xbar[i2])^2 - ϵ  * (+L[2]/2-xbar[i2])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + H[2]/2 - ybar[i2])^2 - ϵ *(H[2]/2 - ybar[i2])^2*(1-W[cur_k]) <= dot_rsqr[i2,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + L[2]/2 - xbar[i2])^2 - ϵ  * (+L[2]/2-xbar[i2])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[2,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - H[2]/2 - ybar[i2])^2 - ϵ *(-H[2]/2 - ybar[i2])^2*(1-W[cur_k]) <= dot_rsqr[i2,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - L[3]/2 - xbar[i3])^2 - ϵ  * (-L[3]/2-xbar[i3])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + H[3]/2 - ybar[i3])^2 - ϵ *(H[3]/2 - ybar[i3])^2*(1-W[cur_k]) <= dot_rsqr[i3,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - L[3]/2 - xbar[i3])^2 - ϵ  * (-L[3]/2-xbar[i3])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - H[3]/2 - ybar[i3])^2 - ϵ *(-H[3]/2 - ybar[i3])^2*(1-W[cur_k]) <= dot_rsqr[i3,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + L[3]/2 - xbar[i3])^2 - ϵ  * (+L[3]/2-xbar[i3])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + H[3]/2 - ybar[i3])^2 - ϵ *(H[3]/2 - ybar[i3])^2*(1-W[cur_k]) <= dot_rsqr[i3,cur_k] )
 				@NLconstraint(s1, ((1-ϵ)*W[cur_k] + ϵ)*(dot_xs[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) + L[3]/2 - xbar[i3])^2 - ϵ  * (+L[3]/2-xbar[i3])^2*(1-W[cur_k]) + ((1-ϵ)*W[cur_k] + ϵ) * (dot_ys[3,cur_k]/((1-ϵ)*W[cur_k] + ϵ) - H[3]/2 - ybar[i3])^2 - ϵ *(-H[3]/2 - ybar[i3])^2*(1-W[cur_k]) <= dot_rsqr[i3,cur_k] )

 			end
 		end
 	end

 	#x=∑ẋ, y=∑ẏ, ∑λ=1
 	@constraint(s1, e1[i in rectangles], sum(dot_x[i,k] for k in max_djc)== x[i])
 	@constraint(s1, e2[i in rectangles], sum(dot_y[i,k] for k in max_djc) == y[i])
 	@constraint(s1, e3[i in rectangles], sum(dot_xs[i,k] for k in max_djc) == xs[i])
 	@constraint(s1, e4[i in rectangles], sum(dot_ys[i,k] for k in max_djc) == ys[i])
 	@constraint(s1, e5[a in areas], sum(dot_rsqr[a,k] for k in max_djc) == rsqr[a])
 	@constraint(s1, e6, sum(W[k] for k in max_djc) == 1)

 	#0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
 	@constraint(s1, b1[i in rectangles,  k in max_djc], dot_x[i,k]<= xub[i] * W[k])
 	@constraint(s1, b2[i in rectangles,  k in max_djc], dot_x[i,k] >= xlb[i] * W[k])
 	@constraint(s1, b3[i in rectangles,  k in max_djc], dot_y[i,k]<= yub[i] * W[k])
 	@constraint(s1, b4[i in rectangles,  k in max_djc], dot_y[i,k] >= ylb[i] * W[k])
 	@constraint(s1, b5[i in rectangles,  k in max_djc], dot_xs[i,k]<= xub[i] * W[k])
 	@constraint(s1, b6[i in rectangles,  k in max_djc], dot_xs[i,k] >= xlb[i] * W[k])
 	@constraint(s1, b7[i in rectangles,  k in max_djc], dot_ys[i,k]<= yub[i] * W[k])
 	@constraint(s1, b8[i in rectangles,  k in max_djc], dot_ys[i,k] >= ylb[i] * W[k])
 	@constraint(s1, b9[a in areas,  k in max_djc], dot_rsqr[a,k] <= 500*W[k])


 	@objective(s1, Min,  prob * sum(rsqr[a] * Carea[a] for a in areas))

 	return s1 
 end