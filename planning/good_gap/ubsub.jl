using JuMP
using Pajarito
using CPLEX 
using Ipopt

function generate_ubsub(; demand=20, qbar=[0 0 0 0 ], ybar=[0 0 0 0], prob=0.25)
    #set up model
    s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0), cont_solver=IpoptSolver(print_level=0)))

    #sets for process
    process = 1:4

    #parameters
    Qu = [2.3 2.8  2 3.2]
    mu = [1 1 1 1]
    α = [80 100 70 110]
    β =0.3* [90 80 100 72]
    ψ = [1 1 1 1]
    γ = 2*[50 50 50 50]
    ϕ=390

    #first stage variables
    @variable(s1, q[i in process]>=0)
    @variable(s1, 0<= y[i in process] <=1)

    #second stage variable
    @variable(s1, R[i in process]>=0)
    @variable(s1, P[i in process]>=0)
    @variable(s1, z[i in process], Bin)
    @variable(s1, δ>=0 )

    # x=x̄
    @constraint(s1, tq[i in process], q[i] == qbar[i])
    @constraint(s1, ty[i in process], y[i] == ybar[i])

    # Ax+g(y) ⩽0, g2(y)⩽0
    @constraint(s1, c1[i in process], z[i] <= y[i])
    @constraint(s1, c2[i in process], P[i]<= Qu[i] * z[i])
    @constraint(s1, c3, sum(P[i] for i in process) == demand - δ)
    @NLconstraint(s1, c4[i in process], -log(1+R[i])* mu[i] + P[i] <= 0)
    @constraint(s1, c5[i in process], P[i] <= q[i])

    #objective
    @objective(s1, Min, prob*(ϕ*δ + sum(α[i]*y[i]+β[i]*q[i]+γ[i]*z[i] + ψ[i] * R[i] for i in process))  )      
    return s1
end