using JuMP
using Pajarito
using CPLEX 
using Ipopt
prob = [0.3 0.4 0.3]
demand = [10 8.5  6]
function generate_fullspace()
    #set up model
    s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0), cont_solver=IpoptSolver(print_level=0)))

    #sets for process
    process = 1:4

    #set for scenarios
    scenarios = 1:3

    #parameters
    Qu = [2.3 2.8  2 3.2]
    mu = [1 1 1 1]
    α = [80 100 70 110]
    β =0.8* [90 80 100 72]
    ψ = [9 9 9 9]
    γ = [50 50 50 50]
    ϕ=290

    #first stage variables
    @variable(s1, q[i in process]>=0)
    @variable(s1, y[i in process], Bin)

    #second stage variable
    @variable(s1, R[i in process, s in scenarios]>=0)
    @variable(s1, P[i in process, s in scenarios]>=0)
    @variable(s1, z[i in process, s in scenarios], Bin)
    @variable(s1, δ[s in scenarios]>=0 )

    #original constraints in the first stage
    @constraint(s1, e1[i in process], q[i]<= y[i] * Qu[i])

    # Ax+g(y) ⩽0, g2(y)⩽0
    @constraint(s1, c1[i in process, s in scenarios], z[i,s] <= y[i])
    @constraint(s1, c2[i in process, s in scenarios], P[i,s]<= Qu[i] * z[i,s])
    @constraint(s1, c3[s in scenarios], sum(P[i,s] for i in process) == demand[s] - δ[s])
    @NLconstraint(s1, c4[i in process, s in scenarios], -log(1+R[i,s])* mu[i] + P[i,s] <= 0)
    @constraint(s1, c5[i in process, s in scenarios], P[i,s] <= q[i])
    #objective
    @objective(s1, Min, sum(prob[s]*(ϕ*δ[s] + sum(α[i]*y[i]+β[i]*q[i]+γ[i]*z[i,s] + ψ[i] * R[i,s] for i in process)) for s in scenarios  ) )     
    return s1
end