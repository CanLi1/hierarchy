using JuMP
using Pajarito
using CPLEX 
using Ipopt
prob = [0.5 0.5]
demand =[1.5 2]
function generate_secondCNF()
    s1 =  Model(solver=PajaritoSolver(rel_gap=0.00001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_THREADS=12), cont_solver=IpoptSolver(print_level=0)))
    scenarios = 1:2
    disjunction = 1:2
    max_djc = 1:2
    ϵ=1e-6 
    @variable(s1, p1>=0)
    @variable(s1, p2>=0)
    @variable(s1, x1, Bin)
    @variable(s1, x2, Bin)

    @variable(s1, y1[s in scenarios]>=0)
    @variable(s1, y2[s in scenarios]>=0)
    @variable(s1, q1[s in scenarios]>=0)
    @variable(s1, q2[s in scenarios] >=0)
    @variable(s1, δ[s in scenarios] >=0)

    @variable(s1, dot_p1[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_p2[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_x1[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_x2[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_y1[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_y2[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_q1[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_q2[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_δ[s in scenarios, d in disjunction, k in max_djc]>=0)
    @variable(s1, dot_λ[s in scenarios, d in disjunction, k in max_djc]>=0)

    #first stage constraints
    @constraint(s1, p1 <= 5 * x1)
    @constraint(s1, p2  <= 6*x2)

    # Ax+g(y) ⩽0, g2(y)⩽0
    @NLconstraint(s1, e1[s in scenarios, d in disjunction, k in max_djc], ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q1[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) -3)^2 - ϵ * 9 *(1-dot_λ[s,d,k]) + ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q2[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) -2)^2 - ϵ * 4 *(1-dot_λ[s,d,k]) <= 17 * dot_λ[s,d,k] -16 * dot_y1[s,d,k])
    @NLconstraint(s1, e2[s in scenarios, d in disjunction, k in max_djc], ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q1[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) -1)^2 - ϵ * 1 *(1-dot_λ[s,d,k]) + ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q2[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) )^2  <= dot_λ[s,d,k] + 16 * dot_y1[s,d,k])
    # @NLconstraint(s1, e3[s in scenarios, d in disjunction, k in max_djc], ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q1[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) -2)^2 - ϵ * 4 *(1-dot_λ[s,d,k]) + ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q2[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) -2)^2 - ϵ * 4 *(1-dot_λ[s,d,k]) <= 17 * dot_λ[s,d,k] -16 * dot_y2[s,d,k])
    # @NLconstraint(s1, e4[s in scenarios, d in disjunction, k in max_djc], ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q1[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ))^2  + ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q2[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) )^2  <= dot_λ[s,d,k] + 16 * dot_y2[s,d,k])
    @NLconstraint(s1, e3[s in scenarios, d in disjunction, k in max_djc], ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q1[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ))^2  + ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q2[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) -1)^2 - ϵ * 1 *(1-dot_λ[s,d,k]) <= 17 * dot_λ[s,d,k] -16 * dot_y2[s,d,k])
    @NLconstraint(s1, e4[s in scenarios, d in disjunction, k in max_djc], ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q1[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ)-4)^2 - ϵ * 16 *(1-dot_λ[s,d,k]) + ((1-ϵ)*dot_λ[s,d,k]+ϵ) * (dot_q2[s,d,k]/((1-ϵ)*dot_λ[s,d,k]+ϵ) -1)^2- ϵ * 1 *(1-dot_λ[s,d,k])  <= dot_λ[s,d,k] + 16 * dot_y2[s,d,k])
    
    @constraint(s1, e5[s in scenarios, d in disjunction, k in max_djc], dot_q1[s,d,k]<= dot_p1[s,d,k])
    @constraint(s1, e6[s in scenarios, d in disjunction, k in max_djc], dot_q2[s,d,k]<= dot_p2[s,d,k])
    @constraint(s1, e7[s in scenarios, d in disjunction, k in max_djc], dot_q1[s,d,k] + dot_q2[s,d,k] + dot_δ[s,d,k] >= demand[s] * dot_λ[s,d,k])
    
    #x=∑ẋ, y=∑ẏ, ∑λ=1
    @constraint(s1, sum1[s in scenarios, d in disjunction], p1 == sum(dot_p1[s,d,k] for k in max_djc))
    @constraint(s1, sum2[s in scenarios, d in disjunction], p2 == sum(dot_p2[s,d,k] for k in max_djc))
    @constraint(s1, sum3[s in scenarios, d in disjunction], x1 == sum(dot_x1[s,d,k] for k in max_djc))
    @constraint(s1, sum4[s in scenarios, d in disjunction], x2 == sum(dot_x2[s,d,k] for k in max_djc))
    @constraint(s1, sum5[s in scenarios, d in disjunction], y1[s] == sum(dot_y1[s,d,k] for k in max_djc))
    @constraint(s1, sum6[s in scenarios, d in disjunction], y2[s] == sum(dot_y2[s,d,k] for k in max_djc))
    @constraint(s1, sum7[s in scenarios, d in disjunction], q1[s] == sum(dot_q1[s,d,k] for k in max_djc))
    @constraint(s1, sum8[s in scenarios, d in disjunction], q2[s] == sum(dot_q2[s,d,k] for k in max_djc))
    @constraint(s1, sum10[s in scenarios, d in disjunction], δ[s] == sum(dot_δ[s,d,k] for k in max_djc))
    @constraint(s1, sum9[s in scenarios, d in disjunction], sum(dot_λ[s,d,k] for k in max_djc) == 1)

    #0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
    @constraint(s1, b1[s in scenarios, d in disjunction, k in max_djc], dot_p1[s, d, k] <= 5*dot_λ[s,d,k])
    @constraint(s1, b2[s in scenarios, d in disjunction, k in max_djc], dot_p2[s, d, k] <= 6*dot_λ[s,d,k])
    @constraint(s1, b3[s in scenarios, d in disjunction, k in max_djc], dot_x1[s, d, k] <= dot_λ[s,d,k])
    @constraint(s1, b4[s in scenarios, d in disjunction, k in max_djc], dot_x2[s, d, k] <= dot_λ[s,d,k])
    @constraint(s1, b5[s in scenarios, d in disjunction, k in max_djc], dot_y1[s, d, k] <= dot_λ[s,d,k])
    @constraint(s1, b6[s in scenarios, d in disjunction, k in max_djc], dot_y2[s, d, k] <= dot_λ[s,d,k])
    @constraint(s1, b7[s in scenarios, d in disjunction, k in max_djc], dot_q1[s, d, k] <= 4*dot_λ[s,d,k])
    @constraint(s1, b8[s in scenarios, d in disjunction, k in max_djc], dot_q2[s, d, k] <= 4*dot_λ[s,d,k])
    @constraint(s1, b9[s in scenarios, d in disjunction, k in max_djc], dot_δ[s, d, k] <= 10*dot_λ[s,d,k])

    #yⱼ=0 or λ
    @constraint(s1, spec1[s in scenarios], dot_y1[s, 1, 1] == dot_λ[s, 1, 1])
    @constraint(s1, spec2[s in scenarios], dot_y1[s, 1, 2] == 0)
    @constraint(s1, spec3[s in scenarios], dot_y2[s, 2, 1] == dot_λ[s, 2, 1])
    @constraint(s1, spec4[s in scenarios], dot_y2[s, 2, 2] == 0)

    @objective(s1, Min, sum(prob[s]* (q1[s]-12*q2[s] + 3*y1[s] -3*y2[s]+ 100 * δ[s]) for s in scenarios ) + 3 * x1 +3 * x2 + p1 + p2)
    return s1
end












