function generate_ubsub(;xbar=zeros(length(rectangles)), ybar=zeros(length(rectangles)), Carea=zeros(areas), prob=0)
 	m = Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0), cont_solver=IpoptSolver(print_level=0)))
 	

 	@variable(m, x[i in rectangles])
 	@variable(m, y[i in rectangles])
 	#second stage var 
 	@variable(m, xs[i in rectangles])
 	@variable(m, ys[i in rectangles])
 	@variable(m, W[i in rectangles, a in areas], Bin)
 	@variable(m, rsqr[a in areas]>=0)

 	#x=x̂
 	@constraint(m, t1[i in rectangles], x[i] == xbar[i])
 	@constraint(m, t2[i in rectangles], y[i] == ybar[i]) 	


 	@constraint(m, c10[i in rectangles, a in areas], (xs[i] - L[i]/2 - xbar[a])^2 + (ys[i] + H[i]/2 - ybar[a])^2 <= rsqr[a] + M2*(1-W[i,a]))
 	@constraint(m, c11[i in rectangles, a in areas], (xs[i] - L[i]/2 - xbar[a])^2 + (ys[i] - H[i]/2 - ybar[a])^2 <= rsqr[a] + M2*(1-W[i,a]))
 	@constraint(m, c12[i in rectangles, a in areas], (xs[i] + L[i]/2 - xbar[a])^2 + (ys[i] + H[i]/2 - ybar[a])^2 <= rsqr[a] + M2*(1-W[i,a]))
 	@constraint(m, c13[i in rectangles, a in areas], (xs[i] + L[i]/2 - xbar[a])^2 + (ys[i] - H[i]/2 - ybar[a])^2 <= rsqr[a] + M2*(1-W[i,a]))
 	@constraint(m, c14[i in rectangles], sum(W[i,a] for a in areas) == 1)
 	@constraint(m, c15[i in rectangles, j in rectangles; i<j], xs[i] - xs[j] == x[i] - x[j])
 	@constraint(m, c16[i in rectangles, j in rectangles; i<j], ys[i] - ys[j] == y[i] - y[j])
 	@objective(m, Min,  prob * sum(rsqr[a] * Carea[a] for a in areas) )
 	return m
 end

