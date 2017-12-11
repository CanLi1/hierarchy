using JuMP
using CPLEX

prob = [0.5 0.5]
demand =[2.5 3]

function generate_master(;mult_p1=[], mult_p2=[],mult_x1=[], mult_x2=[], g=[], iter=[], node1=0, node2=2)
        #set up model
    m1 = Model(solver=CplexSolver(CPX_PARAM_THREADS=1))
    scenarios = 1:2
    @variable(m1, p1>=0)
    @variable(m1, p2>=0)
    @variable(m1, x1, Bin)
    @variable(m1, x2, Bin)

    #node constraint 
    @constraint(m1, p2 >=node1 )
    @constraint(m1, p2<=node2)
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