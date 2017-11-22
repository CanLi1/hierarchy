using JuMP
using CPLEX

prob = [0.25 0.25 0.25 0.25]
Qu = [20 20 20 20]

function generate_master(;mult_q=[], mult_y=[],g=[], iter=[] )

    #set up model
    m1 = Model(solver=CplexSolver())
    scenarios = 1:3

    #sets for process
    process = 1:4

    #continous vars
    @variable(m1, q[i in process] >=0)
    
    #Binary vars
    @variable(m1, y[i in process], Bin)

    #variables for Benders cuts
    @variable(m1, η[s in scenarios] >=0)

    #original constraints in the first stage
    @constraint(m1, e1[i in process], q[i]<= y[i] * Qu[i])

    #Benders cuts 
    @constraint(m1, b1[s in scenarios, it in iter], η[s] >= sum(mult_q[it, s, i] * q[i] + mult_y[it, s, i] * y[i] for i in process ))

    #obj
    @objective(m1, Min, sum(η[s] for s in scenarios))

    return m1 
end