function generate_master(; mult_x=[], mult_y=[], g=[], iter=[], xub = zeros(length(rectangles)), yub = zeros(length(rectangles)), xlb = zeros(length(rectangles)), ylb = zeros(length(rectangles)) )
	m = Model(solver=CplexSolver())

 	@variable(m, delx[i in rectangles, j in rectangles; i<j]>=0)
 	@variable(m, dely[i in rectangles, j in rectangles; i<j]>=0)
 	@variable(m, x[i in rectangles])
 	@variable(m, y[i in rectangles])
 	@variable(m, Z[i in rectangles, j in rectangles, p in pos; i<j], Bin)

 	@variable(m, η[s in scenarios]>=0)
 	@variable(m, firststagecost)

 	@constraint(m, c1[i in rectangles, j in rectangles; i<j], delx[i,j] >= x[i] - x[j])
 	@constraint(m, c2[i in rectangles, j in rectangles; i<j], delx[i,j] >= x[j] - x[i])
 	@constraint(m, c3[i in rectangles, j in rectangles; i<j], dely[i,j] >= y[i] - y[j])
 	@constraint(m, c4[i in rectangles, j in rectangles; i<j], dely[i,j] >= y[j] - y[i])
 	@constraint(m, c5[i in rectangles, j in rectangles; i<j], x[i] + L[i]/2 <= x[j] - L[j] /2 + M*(1 - Z[i,j,1]))
 	@constraint(m, c6[i in rectangles, j in rectangles; i<j], x[j] + L[j]/2 <= x[i] - L[i] /2 + M*(1 - Z[i,j,2]))
 	@constraint(m, c7[i in rectangles, j in rectangles; i<j], y[i] + H[i]/2 <= y[j] - H[j] /2 + M*(1 - Z[i,j,3]))
 	@constraint(m, c8[i in rectangles, j in rectangles; i<j], y[j] + H[j]/2 <= y[i] - H[i] /2 + M*(1 - Z[i,j,4]))
 	@constraint(m, c9[i in rectangles, j in rectangles; i<j], sum(Z[i,j,p] for p in pos) == 1)

 	#bound on first stage variables 
 	@constraint(m, bx1[i in rectangles], x[i]>=xlb[i])
 	@constraint(m, bx2[i in rectangles], x[i]<=xub[i])
 	@constraint(m, by1[i in rectangles], y[i]>=ylb[i])
 	@constraint(m, by2[i in rectangles], y[i]<=yub[i])
 	#Benders cuts
 	@constraint(m, b1[s in scenarios, it in iter], η[s]>= g[it][s] + sum(mult_x[it][s][i] * x[i] + mult_y[it][s][i]*y[i] for i in rectangles))

 	#first stage cost 
 	@constraint(m, firststagecost== sum(C[i,j] * (dely[i,j] + delx[i,j]) for i in rectangles for j in rectangles if i<j))
 	@objective(m, Min, sum(η[s] for s in scenarios) + sum(C[i,j] * (dely[i,j] + delx[i,j]) for i in rectangles for j in rectangles if i<j) )

 	return m
 end