using Mosek
function generate_mosek_sub(; p1bar=0, p2bar=0, x1bar=0, x2bar=0, prob=0.5, demand=2, ub=[4 2], lb=[0 0])
    # s1 =  Model(solver=MosekSolver(MSK_IPAR_NUM_THREADS=1))
    s1 =  Model(solver=MosekSolver(MSK_DPAR_OPTIMIZER_MAX_TIME=50.0, MSK_IPAR_NUM_THREADS=1, MSK_IPAR_INTPNT_MAX_ITERATIONS=1000000, MSK_DPAR_INTPNT_TOL_REL_GAP=1e-5, MSK_DPAR_INTPNT_NL_TOL_REL_GAP=1e-5))
    # s1 = Model(solver=IpoptSolver())
    # s1 = Model(solver=KnitroSolver())
    max_djc = 1:4
    ϵ=1e-6 
    @variable(s1, p1>=0)
    @variable(s1, p2>=0)
    @variable(s1, x1)
    @variable(s1, x2)

    @variable(s1, y1>=0)
    @variable(s1, y2>=0)
    @variable(s1, q1>=0)
    @variable(s1, q2 >=0)
    @variable(s1, δ>=0)

    @variable(s1, dot_p1[k in max_djc]>=0)
    @variable(s1, dot_p2[k in max_djc]>=0)
    @variable(s1, dot_x1[k in max_djc]>=0)
    @variable(s1, dot_x2[k in max_djc]>=0)
    @variable(s1, dot_y1[k in max_djc]>=0)
    @variable(s1, dot_y2[k in max_djc]>=0)
    @variable(s1, dot_q1[k in max_djc]>=0)
    @variable(s1, dot_q2[k in max_djc]>=0)
    @variable(s1, dot_δ[k in max_djc]>=0)
    @variable(s1, dot_λ[k in max_djc]>=0)

    @constraint(s1, t1, p1== p1bar)
    @constraint(s1, t2, p2== p2bar)
    @constraint(s1, t3, x1 == x1bar)
    @constraint(s1, t4, x2 == x2bar)
 # Ax+g(y) ⩽0, g2(y)⩽0
    @NLconstraint(s1, e1[k in max_djc], ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q1[k]/((1-ϵ)*dot_λ[k]+ϵ) -3)^2 - ϵ * 9 *(1-dot_λ[k]) + ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q2[k]/((1-ϵ)*dot_λ[k]+ϵ) -2)^2 - ϵ * 4 *(1-dot_λ[k]) <= 17 * dot_λ[k] -16 * dot_y1[k])
    @NLconstraint(s1, e2[k in max_djc], ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q1[k]/((1-ϵ)*dot_λ[k]+ϵ) -1)^2 - ϵ * 1 *(1-dot_λ[k]) + ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q2[k]/((1-ϵ)*dot_λ[k]+ϵ) )^2  <= dot_λ[k] + 16 * dot_y1[k])
    @NLconstraint(s1, e3[k in max_djc], ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q1[k]/((1-ϵ)*dot_λ[k]+ϵ))^2  + ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q2[k]/((1-ϵ)*dot_λ[k]+ϵ) -1)^2 - ϵ * 1 *(1-dot_λ[k]) <= 17 * dot_λ[k] -16 * dot_y2[k])
    @NLconstraint(s1, e4[k in max_djc], ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q1[k]/((1-ϵ)*dot_λ[k]+ϵ)-4)^2 - ϵ * 16 *(1-dot_λ[k]) + ((1-ϵ)*dot_λ[k]+ϵ) * (dot_q2[k]/((1-ϵ)*dot_λ[k]+ϵ) -1)^2- ϵ * 1 *(1-dot_λ[k])  <= dot_λ[k] + 16 * dot_y2[k])
    
    @constraint(s1, e5[k in max_djc], dot_q1[k]<= dot_p1[k])
    @constraint(s1, e6[k in max_djc], dot_q2[k]<= dot_p2[k])
    @constraint(s1, e7[k in max_djc], dot_q1[k] + dot_q2[k] + dot_δ[k] >= demand * dot_λ[k])
    
    #x=∑ẋ, y=∑ẏ, ∑λ=1
    @constraint(s1, sum1, p1 == sum(dot_p1[k] for k in max_djc))
    @constraint(s1, sum2, p2 == sum(dot_p2[k] for k in max_djc))
    @constraint(s1, sum3, x1 == sum(dot_x1[k] for k in max_djc))
    @constraint(s1, sum4, x2 == sum(dot_x2[k] for k in max_djc))
    @constraint(s1, sum5, y1 == sum(dot_y1[k] for k in max_djc))
    @constraint(s1, sum6, y2 == sum(dot_y2[k] for k in max_djc))
    @constraint(s1, sum7, q1 == sum(dot_q1[k] for k in max_djc))
    @constraint(s1, sum8, q2 == sum(dot_q2[k] for k in max_djc))
    @constraint(s1, sum10, δ == sum(dot_δ[k] for k in max_djc))
    @constraint(s1, sum9, sum(dot_λ[k] for k in max_djc) == 1)

    #0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
    @constraint(s1, b1[k in max_djc], dot_p1[k] <= ub[1] * dot_λ[k])
    @constraint(s1, b2[k in max_djc], dot_p2[k] <= ub[2]*dot_λ[k])
    @constraint(s1, bb1[k in max_djc], dot_p1[k] >= lb[1]*dot_λ[k])
    @constraint(s1, bb2[k in max_djc], dot_p2[k] >= lb[2]*dot_λ[k])
    @constraint(s1, b3[k in max_djc], dot_x1[k] <= dot_λ[k])
    @constraint(s1, b4[k in max_djc], dot_x2[k] <= dot_λ[k])
    @constraint(s1, b5[k in max_djc], dot_y1[k] <= dot_λ[k])
    @constraint(s1, b6[k in max_djc], dot_y2[k] <= dot_λ[k])
    @constraint(s1, b7[k in max_djc], dot_q1[k] <= 4*dot_λ[k])
    @constraint(s1, b8[k in max_djc], dot_q2[k] <= 2*dot_λ[k])
    @constraint(s1, b9[k in max_djc], dot_δ[k] <= 10*dot_λ[k])

    #yⱼ=0 or λ
    @constraint(s1, spec1, dot_y1[ 1] == dot_λ[ 1])
    @constraint(s1, spec2, dot_y1[ 2] == 0)
    @constraint(s1, spec3, dot_y1[ 3] == dot_λ[ 3])
    @constraint(s1, spec4, dot_y1[ 4] == 0)
    @constraint(s1, spec5, dot_y2[ 1] == dot_λ[ 1])
    @constraint(s1, spec6, dot_y2[ 2] == 0)
    @constraint(s1, spec7, dot_y2[ 3] == 0)
    @constraint(s1, spec8, dot_y2[ 4] == dot_λ[4])

    @objective(s1, Min, prob* (q1-12*q2 + 3*y1 -3*y2 + 3 * x1 +3 * x2 + p1 + p2 + 100 * δ))
    return s1
end












