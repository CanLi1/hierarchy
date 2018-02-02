 using JuMP
 using Pajarito
 using CPLEX 
 using Ipopt
 using Mosek
 using BARON
include("input.jl")
 function generate_fullspace()
 	m = Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(), cont_solver=IpoptSolver()))
 	# m = Model(solver=BaronSolver())
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

 	djc1 = 1:4
 	djc2 = 1:2
 	ϵ = 1e-6

 	@variable(m, dot_x[i in rectangles, ii in rectangles, jj in rectangles, k in djc1; ii<jj])
 	@variable(m, dot_y[i in rectangles, ii in rectangles, jj in rectangles, k in djc1; ii<jj])
 	@variable(m, dot_xs[i in rectangles, s in scenarios, a in areas])
 	@variable(m, dot_ys[i in rectangles, s in scenarios, a in areas])
 	@variable(m, dot_rsqr[a in areas, s in scenarios, i in rectangles, aa in areas ]>=0)

 	@constraint(m, c5[i in rectangles, j in rectangles; i<j], dot_x[i, i, j, 1] + L[i]/2 * Z[i,j, 1] <= dot_x[j, i, j, 1] - L[j] /2 * Z[i,j,1])
 	@constraint(m, c6[i in rectangles, j in rectangles; i<j], dot_x[j, i, j, 2] + L[j]/2 * Z[i,j,2] <= dot_x[i, i,j,2] - L[i] /2 *Z[i,j,2])
 	@constraint(m, c7[i in rectangles, j in rectangles; i<j], dot_y[i,i,j,3] + H[i]/2 *Z[i,j,3] <= dot_y[j,i,j,3] - H[j] /2 *Z[i,j,3])
 	@constraint(m, c8[i in rectangles, j in rectangles; i<j], dot_y[j,i,j,4] + H[j]/2*Z[i,j,4] <= dot_y[i,i,j,4] - H[i] /2 *Z[i,j,4])
 	@constraint(m, c9[i in rectangles, j in rectangles; i<j], sum(Z[i,j,p] for p in pos) == 1)

 	@constraint(m, e1[i in rectangles, ii in rectangles, jj in rectangles; ii<jj], x[i] == sum(dot_x[i, ii, jj, p] for p in pos))
 	@constraint(m, e2[i in rectangles, ii in rectangles, jj in rectangles; ii<jj], y[i] == sum(dot_y[i, ii, jj, p] for p in pos))
 	@constraint(m, b1[i in rectangles, ii in rectangles, jj in rectangles, p in pos;ii<jj], dot_x[i, ii, jj, p]<= Z[ii,jj,p] * 10)
 	@constraint(m, b2[i in rectangles, ii in rectangles, jj in rectangles, p in pos;ii<jj], dot_x[i, ii, jj, p]>= Z[ii,jj,p] * (-10))
 	@constraint(m, b3[i in rectangles, ii in rectangles, jj in rectangles, p in pos;ii<jj], dot_y[i, ii, jj, p]<= Z[ii,jj,p] * 40)
 	@constraint(m, b4[i in rectangles, ii in rectangles, jj in rectangles, p in pos;ii<jj], dot_y[i, ii, jj, p]>= Z[ii,jj,p] * (-5)) 


 	@constraint(m, c1[i in rectangles, j in rectangles; i<j], delx[i,j] >= x[i] - x[j])
 	@constraint(m, c2[i in rectangles, j in rectangles; i<j], delx[i,j] >= x[j] - x[i])
 	@constraint(m, c3[i in rectangles, j in rectangles; i<j], dely[i,j] >= y[i] - y[j])
 	@constraint(m, c4[i in rectangles, j in rectangles; i<j], dely[i,j] >= y[j] - y[i])

 	@constraint(m, c15[i in rectangles, j in rectangles, s in scenarios; i<j], xs[i,s] - xs[j,s] == x[i] - x[j])
 	@constraint(m, c16[i in rectangles, j in rectangles, s in scenarios; i<j], ys[i,s] - ys[j,s] == y[i] - y[j])



 	@NLconstraint(m, c10[i in rectangles, a in areas, s in scenarios], ((1-ϵ)*W[i,a,s] + ϵ)*(dot_xs[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) - L[i]/2 - xbar[a])^2 - ϵ  * (-L[i]/2-xbar[a])^2*(1-W[i,a,s]) + ((1-ϵ)*W[i,a,s] + ϵ) * (dot_ys[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) + H[i]/2 - ybar[a])^2 - ϵ *(H[i]/2 - ybar[a])^2*(1-W[i,a,s]) <= dot_rsqr[a,s,i,a] )
 	@NLconstraint(m, c11[i in rectangles, a in areas, s in scenarios], ((1-ϵ)*W[i,a,s] + ϵ)*(dot_xs[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) - L[i]/2 - xbar[a])^2 - ϵ  * (-L[i]/2-xbar[a])^2*(1-W[i,a,s]) + ((1-ϵ)*W[i,a,s] + ϵ) * (dot_ys[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) - H[i]/2 - ybar[a])^2 - ϵ *(-H[i]/2 - ybar[a])^2*(1-W[i,a,s]) <= dot_rsqr[a,s,i,a] )
 	@NLconstraint(m, c12[i in rectangles, a in areas, s in scenarios], ((1-ϵ)*W[i,a,s] + ϵ)*(dot_xs[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) + L[i]/2 - xbar[a])^2 - ϵ  * (+L[i]/2-xbar[a])^2*(1-W[i,a,s]) + ((1-ϵ)*W[i,a,s] + ϵ) * (dot_ys[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) + H[i]/2 - ybar[a])^2 - ϵ *(H[i]/2 - ybar[a])^2*(1-W[i,a,s]) <= dot_rsqr[a,s,i,a] )
 	@NLconstraint(m, c13[i in rectangles, a in areas, s in scenarios], ((1-ϵ)*W[i,a,s] + ϵ)*(dot_xs[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) + L[i]/2 - xbar[a])^2 - ϵ  * (+L[i]/2-xbar[a])^2*(1-W[i,a,s]) + ((1-ϵ)*W[i,a,s] + ϵ) * (dot_ys[i,s,a]/((1-ϵ)*W[i,a,s] + ϵ) - H[i]/2 - ybar[a])^2 - ϵ *(-H[i]/2 - ybar[a])^2*(1-W[i,a,s]) <= dot_rsqr[a,s,i,a] )
 	@constraint(m, e3[i in rectangles, s in scenarios], sum(dot_xs[i,s,a] for a in areas) == xs[i,s])
 	@constraint(m, e4[i in rectangles, s in scenarios], sum(dot_ys[i,s,a] for a in areas) == ys[i,s])
 	@constraint(m, e5[a in areas, s in scenarios, i in rectangles], rsqr[a,s] == sum(dot_rsqr[a,s,i,aa] for aa in areas ))
 	@constraint(m, e6[i in rectangles, s in scenarios], sum(W[i,a,s] for a in areas) == 1)
 	@constraint(m, b5[i in rectangles, s in scenarios, a in areas], dot_xs[i,s,a] <= 10*W[i,a,s])
 	@constraint(m, b6[i in rectangles, s in scenarios, a in areas], dot_xs[i,s,a] >= (-10)*W[i,a,s])
 	@constraint(m, b7[i in rectangles, s in scenarios, a in areas], dot_ys[i,s,a] <= 40*W[i,a,s])
 	@constraint(m, b8[i in rectangles, s in scenarios, a in areas], dot_ys[i,s,a] >= (-5)*W[i,a,s])
 	@constraint(m, b9[a in areas, s in scenarios, i in rectangles, aa in areas], dot_rsqr[a,s,i,aa] <= 500 * W[i,aa,s])


 	@objective(m, Min, sum(C[i,j] * (dely[i,j] + delx[i,j]) for i in rectangles for j in rectangles if i<j) + sum(prob[s] * sum(rsqr[a,s] * Carea[a,s] for a in areas) for s in scenarios))
 	# @objective(m, Min, sum(C[i,j] * (dely[i,j] + delx[i,j]) for i in rectangles for j in rectangles if i<j) )
 	return m
 end

m = generate_fullspace()

 solve(m)
 println(getvalue(getindex(m, :x)))
 println(getvalue(getindex(m, :y)))
 println(getvalue(getindex(m, :Z)))
 # println(getvalue(getindex(m, :dot_x)))
 # println(getvalue(getindex(m, :dot_y)))
 println(getvalue(getindex(m, :W)))
 println(getvalue(getindex(m, :rsqr)))
 # println(getvalue(getindex(m, :xs)))
 # println(getvalue(getindex(m, :ys)))
 #  println(getvalue(getindex(m, :dot_xs)))
 # println(getvalue(getindex(m, :dot_ys)))
 # println(getvalue(getindex(m, :dot_rsqr)))