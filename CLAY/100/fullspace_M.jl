 using JuMP
 using Pajarito
 using CPLEX 
 using Ipopt
include("input.jl")
 function generate_fullspace()
 	m = Model(solver=PajaritoSolver(timeout=10000, rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0), cont_solver=IpoptSolver(print_level=0)))
 	
 	@variable(m, delx[i in rectangles, j in rectangles; i<j]>=0)
 	@variable(m, dely[i in rectangles, j in rectangles; i<j]>=0)
 	@variable(m, x[i in rectangles])
 	@variable(m, y[i in rectangles])
 	@variable(m, Z[i in rectangles, j in rectangles, p in pos; i<j], Bin)

 	#second stage var 
 	@variable(m, xs[i in rectangles, s in scenarios])
 	@variable(m, ys[i in rectangles, s in scenarios])
 	@variable(m, W[i in rectangles, a in areas, s in scenarios], Bin)
 	@variable(m, rsqr[a in areas, s in scenarios]>=0)

 	@constraint(m, c1[i in rectangles, j in rectangles; i<j], delx[i,j] >= x[i] - x[j])
 	@constraint(m, c2[i in rectangles, j in rectangles; i<j], delx[i,j] >= x[j] - x[i])
 	@constraint(m, c3[i in rectangles, j in rectangles; i<j], dely[i,j] >= y[i] - y[j])
 	@constraint(m, c4[i in rectangles, j in rectangles; i<j], dely[i,j] >= y[j] - y[i])
 	@constraint(m, c5[i in rectangles, j in rectangles; i<j], x[i] + L[i]/2 <= x[j] - L[j] /2 + M*(1 - Z[i,j,1]))
 	@constraint(m, c6[i in rectangles, j in rectangles; i<j], x[j] + L[j]/2 <= x[i] - L[i] /2 + M*(1 - Z[i,j,2]))
 	@constraint(m, c7[i in rectangles, j in rectangles; i<j], y[i] + H[i]/2 <= y[j] - H[j] /2 + M*(1 - Z[i,j,3]))
 	@constraint(m, c8[i in rectangles, j in rectangles; i<j], y[j] + H[j]/2 <= y[i] - H[i] /2 + M*(1 - Z[i,j,4]))
 	@constraint(m, c9[i in rectangles, j in rectangles; i<j], sum(Z[i,j,p] for p in pos) == 1)

 	@constraint(m, c10[i in rectangles, a in areas, s in scenarios], (xs[i,s] - L[i]/2 - xbar[a])^2 + (ys[i,s] + H[i]/2 - ybar[a])^2 <= rsqr[a,s] + M2*(1-W[i,a,s]))
 	@constraint(m, c11[i in rectangles, a in areas, s in scenarios], (xs[i,s] - L[i]/2 - xbar[a])^2 + (ys[i,s] - H[i]/2 - ybar[a])^2 <= rsqr[a,s] + M2*(1-W[i,a,s]))
 	@constraint(m, c12[i in rectangles, a in areas, s in scenarios], (xs[i,s] + L[i]/2 - xbar[a])^2 + (ys[i,s] + H[i]/2 - ybar[a])^2 <= rsqr[a,s] + M2*(1-W[i,a,s]))
 	@constraint(m, c13[i in rectangles, a in areas, s in scenarios], (xs[i,s] + L[i]/2 - xbar[a])^2 + (ys[i,s] - H[i]/2 - ybar[a])^2 <= rsqr[a,s] + M2*(1-W[i,a,s]))
 	@constraint(m, c14[i in rectangles, s in scenarios], sum(W[i,a,s] for a in areas) == 1)
 	@constraint(m, c15[i in rectangles, j in rectangles, s in scenarios; i<j], xs[i,s] - xs[j,s] == x[i] - x[j])
 	@constraint(m, c16[i in rectangles, j in rectangles, s in scenarios; i<j], ys[i,s] - ys[j,s] == y[i] - y[j])
 	@objective(m, Min, sum(C[i,j] * (dely[i,j] + delx[i,j]) for i in rectangles for j in rectangles if i<j) + sum(prob[s] * sum(rsqr[a,s] * Carea[a,s] for a in areas) for s in scenarios))
 	return m
 end

 m = generate_fullspace()

 solve(m)
 println(getvalue(getindex(m, :x)))
 println(getvalue(getindex(m, :y)))
 println(getvalue(getindex(m, :Z)))
 println(getvalue(getindex(m, :W)))
 println(getvalue(getindex(m, :rsqr)))
