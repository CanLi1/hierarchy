using JuMP
using Mosek


function generate_sub(;djc=[[1,2,3,4]], demand=20.0, qbar=[0.0 0.0 0.0 0.0 ], ybar=[0.0 0.0 0.0 0.0], prob=0.25)
    #set up model
    s1 =  Model(solver=MosekSolver(MSK_IPAR_NUM_THREADS=0))
    # s1 = Model(solver=CplexSolver())
    ϵ=1e-5
    #sets for number of disjunctions
    disjunction=1:length(djc)

    #sets for process
    process = 1:4
    
    #maximum number of disjuncts 
    max_djc = 1:16

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
    @variable(s1, 0<= y[i in process] <=1)

    #second stage variable
    @variable(s1, R[i in process]>=0)
    @variable(s1, P[i in process]>=0)
    @variable(s1, 0<=z[i in process]<=1)
    @variable(s1, δ>=0 )

    #variable for disjunction 
    @variable(s1, dot_q[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_y[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_R[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_P[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_z[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_δ[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)
    @variable(s1, dot_λ[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))]>=0)

    # x=x̄
    @constraint(s1, tq[i in process], q[i] == qbar[i])
    @constraint(s1, ty[i in process], y[i] == ybar[i])

    # Ax+g(y) ⩽0, g2(y)⩽0
    @constraint(s1, c1[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_z[i,d,k] <= dot_y[i,d,k])
    @constraint(s1, c2[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_P[i,d,k]<= Qu[i] * dot_z[i,d,k])
    @constraint(s1, c3[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], sum(dot_P[i,d,k] for i in process) == demand*dot_λ[d,k] - dot_δ[d,k])
    @NLconstraint(s1, c4[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], -( (1-ϵ)* dot_λ[d,k]  + ϵ)*log(1+dot_R[i,d,k]/((1-ϵ)*dot_λ[d,k] + ϵ))* mu[i] + dot_P[i,d,k] <= 0)
    @constraint(s1, c5[i in process,  d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_P[i,d,k] <= dot_q[i,d,k])
    
    #x=∑ẋ, y=∑ẏ, ∑λ=1
    @constraint(s1, e1[i in process, d in disjunction], q[i] == sum(dot_q[i,d,k] for k in max_djc if k<= 2^(length(djc[d])) ))
    @constraint(s1, e2[i in process, d in disjunction], y[i] == sum(dot_y[i,d,k] for k in max_djc if k<= 2^(length(djc[d]))))
    @constraint(s1, e3[i in process, d in disjunction], R[i] == sum(dot_R[i,d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e4[i in process, d in disjunction], P[i] == sum(dot_P[i,d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e5[i in process, d in disjunction], z[i] == sum(dot_z[i,d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e6[d in disjunction], δ == sum(dot_δ[d,k] for k in max_djc if k <=2^(length(djc[d]))))
    @constraint(s1, e7[d in disjunction], sum(dot_λ[d,k] for k in max_djc if k <=2^(length(djc[d]))) == 1)

    #0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
    @constraint(s1, b1[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_q[i,d,k] <= Qu[i] * dot_λ[d,k])
    @constraint(s1, b2[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_y[i,d,k] <= dot_λ[d,k])
    @constraint(s1, b3[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_P[i,d,k] <= Qu[i] * dot_λ[d,k])
    @constraint(s1, b4[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_R[i,d,k] <= (exp(Qu[i])-1)* dot_λ[d,k])
    @constraint(s1, b5[i in process, d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_z[i,d,k] <= dot_λ[d,k])
    @constraint(s1, b6[d in disjunction, k in max_djc; k<= 2^(length(djc[d]))], dot_δ[d,k] <= demand * dot_λ[d,k])
    
    #yⱼ=0 or 1
    index_d = 1
    for item in djc
        for index_k in 1:(2^length(item))
            seq_i = 1
            for index_i in item
                if (floor(Int, (index_k-1) / 2^(seq_i-1) ) % 2) == 1
                    @constraint(s1, dot_z[index_i, index_d, index_k] == dot_λ[index_d, index_k])
                else
                    @constraint(s1, dot_z[index_i, index_d, index_k] == 0)
                end
                seq_i = seq_i + 1
            end 
        end
        index_d = index_d + 1
    end

    #objective
    @objective(s1, Min, prob*(ϕ*δ + sum(α[i]*y[i]+β[i]*q[i]+γ[i]*z[i] + ψ[i] * R[i] for i in process))  )      
    return s1
end 