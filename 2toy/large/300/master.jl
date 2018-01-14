using CPLEX
function generate_master(;mult_p1=[], mult_p2=[],mult_x1=[], mult_x2=[], g=[], iter=[], lb=[0 0], ub=[4 2])
        #set up model
    m1 = Model(solver=CplexSolver())
    @variable(m1, lb[1] <=p1<=ub[1])
    @variable(m1, lb[2] <= p2<=ub[2])
    @variable(m1, x1, Bin)
    @variable(m1, x2, Bin)

    #node constraint 
    @variable(m1, η[s in scenarios]>=-100)
    #first stage constraints
    @constraint(m1, p1 <= 4 * x1)
    @constraint(m1, p2  <= 2*x2)  
    @constraint(m1, init[s in scenarios], η[s] >= -100 +0.5*( 3 * x1 +3 * x2 + p1 + p2))
    #Benders cuts 
    @constraint(m1, b1[s in scenarios, i in iter], η[s] >= g[s,i] + mult_p1[s,i] * p1 + mult_p2[s,i]*p2 +mult_x1[s,i] * x1 +mult_x2[s,i] * x2  )
    #objective
    @objective(m1, Min, sum(η[s] for s in scenarios))
    return m1
end