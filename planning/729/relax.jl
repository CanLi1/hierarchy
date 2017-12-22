using JuMP
using Pajarito
using CPLEX 
using Mosek
# using Ipopt
# using KNITRO
baseprob = [0.3 0.4 0.3]
basedemand = [10 8.5  6]
basemu = [0.9 1 1.1]
baseϕ = [360 390 420]
prob = zeros(Float64, 243)
demand = zeros(Float64, 243)
mu=zeros(Float64, 4, 243)
ϕ =zeros(Float64, 243)
for i2 in 1:3
    for i3 in 1:3
        for i4 in 1:3
            for i5 in 1:3
                for i6 in 1:3
                    temp_scenario = 81 *(i2-1) + 27*(i3-1) + 9*(i4-1) + 3*(i5-1) + i6
                    demand[temp_scenario] = basedemand[i2]
                    mu[1, temp_scenario] = basemu[i3]
                    mu[2, temp_scenario] = basemu[i4]
                    mu[3, temp_scenario] = basemu[i5]
                    mu[4, temp_scenario] = basemu[i5]
                    ϕ[temp_scenario] = baseϕ[i6]
                    prob[temp_scenario] =  baseprob[i2] * baseprob[i3] * baseprob[i4] * baseprob[i5] * baseprob[i6]
                end
            end
        end
    end
end

function generate_fullspace()
    #set up model
    s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(), cont_solver=MosekSolver(MSK_IPAR_NUM_THREADS=0)))
    # s1 = Model(solver=IpoptSolver())

    #sets for process
    process = 1:4

    #set for scenarios
    scenarios = 1:243

    #parameters
    Qu = [2.3 2.8  2 3.2]
    α = [80 100 70 110]
    β =0.3* [90 80 100 72]
    ψ = [1 1 1 1]
    γ = 2*[50 50 50 50]


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
    @NLconstraint(s1, c4[i in process, s in scenarios], -log(1+R[i,s])* mu[i,s] + P[i,s] <= 0)
    @constraint(s1, c5[i in process, s in scenarios], P[i,s] <= q[i])
    #objective
    @objective(s1, Min, sum(prob[s]*(ϕ[s]*δ[s] + sum(α[i]*y[i]+β[i]*q[i]+γ[i]*z[i,s] + ψ[i] * R[i,s] for i in process)) for s in scenarios  ) )     
    return s1
end