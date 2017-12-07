using JuMP
using Pajarito
using CPLEX 
using Ipopt
prob = [0.3 0.4 0.3]
demand = [10 8.5  6]
function generate_fullspaceCH(;ϕ=200, β=[90 80 100 72], ψ=[1 1 1 1], γ=[50 50 50 50])
    #set up model
    s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_THREADS=12), cont_solver=IpoptSolver(print_level=0)))
    djc=[[1],[2],[3],[4]]
    ϵ=1e-5    
    #sets for number of disjunctions
    disjunction=1:length(djc)
    #maximum number of disjuncts 
    max_djc = 1:16
    #sets for process
    process = 1:4

    #set for scenarios
    scenarios = 1:3

    #parameters
    Qu = [2.3 2.8  2 3.2]
    mu = [1 1 1 1]
    α = [80 100 70 110]



    #first stage variables
    @variable(s1, q[i in process]>=0)
    @variable(s1, y[i in process], Bin)

    #second stage variable
    @variable(s1, R[ s in scenarios, i in process]>=0)
    @variable(s1, P[ s in scenarios, i in process ]>=0)
    @variable(s1, z[ s in scenarios, i in process])
    @variable(s1, δ[s in scenarios]>=0 )

    #variable for disjunction 
    @variable(s1, dot_q[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_y[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_R[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_P[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_z[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_δ[s in scenarios, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_λ[s in scenarios, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)


    #original constraints in the first stage
    @constraint(s1, f1[i in process], q[i]<= y[i] * Qu[i])

# Ax+g(y) ⩽0, g2(y)⩽0
    @constraint(s1, c1[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_z[s,i,d,k] <= dot_y[s,i,d,k])
    @constraint(s1, c2[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_P[s,i,d,k]<= Qu[i] * dot_z[s,i,d,k])
    @constraint(s1, c3[s in scenarios, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], sum(dot_P[s,i,d,k] for i in process) == demand[s]*dot_λ[s,d,k] - dot_δ[s,d,k])
    @NLconstraint(s1, c4[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], -( (1-ϵ)* dot_λ[s,d,k]  + ϵ)*log(1+dot_R[s,i,d,k]/((1-ϵ)*dot_λ[s,d,k] + ϵ))* mu[i] + dot_P[s,i,d,k] <= 0)
    @constraint(s1, c5[s in scenarios, i in process,  d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_P[s,i,d,k] <= dot_q[s,i,d,k])
    
    #x=∑ẋ, y=∑ẏ, ∑λ=1
    @constraint(s1, e1[s in scenarios, i in process, d in disjunction], q[i] == sum(dot_q[s, i,d,k] for k in max_djc if k<= 2^(length(djc[d])) ))
    @constraint(s1, e2[s in scenarios, i in process, d in disjunction], y[i] == sum(dot_y[s, i,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e3[s in scenarios, i in process, d in disjunction], R[s,i] == sum(dot_R[s,i,d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e4[s in scenarios, i in process, d in disjunction], P[s,i] == sum(dot_P[s, i,d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e5[s in scenarios, i in process, d in disjunction], z[s,i] == sum(dot_z[s,i,d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e6[s in scenarios, d in disjunction], δ[s] == sum(dot_δ[s,d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e7[s in scenarios, d in disjunction], sum(dot_λ[s,d,k] for k in max_djc if k <=2^(length(djc[d]))) == 1)

    #0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
    @constraint(s1, b1[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_q[s,i,d,k] <= Qu[i] * dot_λ[s,d,k])
    @constraint(s1, b2[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_y[s,i,d,k] <= dot_λ[s,d,k])
    @constraint(s1, b3[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_P[s,i,d,k] <= Qu[i] * dot_λ[s,d,k])
    @constraint(s1, b4[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_R[s,i,d,k] <= (exp(Qu[i])-1)* dot_λ[s,d,k])
    @constraint(s1, b5[s in scenarios, i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_z[s,i,d,k] <= dot_λ[s,d,k])
    @constraint(s1, b6[s in scenarios, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_δ[s,d,k] <= demand[s] * dot_λ[s,d,k])
    
    #yⱼ=0 or 1
    
    for s in scenarios
        index_d = 1
        for item in djc
            for index_k in 1:(2^length(item))
                seq_i = 1
                for index_i in item
                    if (floor(Int, (index_k-1) / 2^(seq_i-1) ) % 2) == 1
                        @constraint(s1, dot_z[s,index_i, index_d, index_k] == dot_λ[s,index_d, index_k])
                    else
                        @constraint(s1,  dot_z[s, index_i, index_d, index_k] == 0)
                    end
                    seq_i = seq_i + 1
                end 
            end
            index_d = index_d + 1
        end
    end
    #objective
    @objective(s1, Min, sum(prob[s]*(ϕ*δ[s] + sum(α[i]*y[i]+β[i]*q[i]+γ[i]*z[s,i] + ψ[i] * R[s,i] for i in process)) for s in scenarios  ) )     
    return s1
end