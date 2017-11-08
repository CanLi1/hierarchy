using JuMP
using CPLEX
prob = [0.5 0.5]
function generate_master(;mult_u=[], mult_v=[],g=[], iter=[] )

    #set up model
    m1 = Model(solver=CplexSolver())
    scenarios = 1:2 
    #continous vars
    @variable(m1, qu >=0)
    @variable(m1, qv>=0)
    @variable(m1, first_stage_cost)
    #binary variable
    @variable(m1, y1, Bin)
    @variable(m1, y2, Bin)

    #variables for Benders cuts
    @variable(m1, yita[s in scenarios] >= -500)

    #original constraint in the first stage
    @constraint(m1, e1, qu <= 2 * y1)
    @constraint(m1, e2, qv <= 3*y2)

    @constraint(m1, first_stage_cost == sum(prob[s]*(10 * y1 +14 * y2+2* (qu+qv) ) for s in scenarios))
    #Benders cuts
    @constraint(m1, b1[s in scenarios, i in iter], yita[s]>= prob[s]*(10 * y1 +14 * y2+2* (qu+qv) )+ mult_u[s,i] * qu + mult_v[s,i] * qv + g[s,i] )

    #objective
    @objective(m1, Min, sum(yita[s] for s in scenarios))

    return m1 
end