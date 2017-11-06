using JuMP
using Ipopt

function generate_sub(;qubar=0, qvbar=0, probs=0.5, ds=20)
    #set up model
    s1 =  Model(solver=IpoptSolver(print_level=0))
    #sets for disjunct 
    lz1 = 1:2
    lz2 = 1:2
    ϵ=1e-5
    # @NLparameter(s1,  qub == qubar)
    # @NLparameter(s1,  qvb == qvbar)

    #continous vars
    @variable(s1, u>=0)
    @variable(s1, v>=0)
    @variable(s1, up>=0)
    @variable(s1, vp>=0)
    @variable(s1, qu)
    @variable(s1, qv)

    #relaxed binary vars
    @variable(s1, 0<=z1<=1)
    @variable(s1, 0<=z2<=1)

    #variables for each disjunct
    @variable(s1, dot_u[i in lz1,j in lz2]>=0)
    @variable(s1, dot_v[i in lz1,j in lz2]>=0)
    @variable(s1, dot_up[i in lz1,j in lz2]>=0)
    @variable(s1, dot_vp[i in lz1,j in lz2]>=0)
    @variable(s1, dot_qu[i in lz1,j in lz2]>=0)
    @variable(s1, dot_qv[i in lz1,j in lz2]>=0)

    #relaxed binary vars
    @variable(s1, 0<=dot_z1[i in lz1,j in lz2]<=1)
    @variable(s1, 0<=dot_z2[i in lz1,j in lz2]<=1)

    #variable that determine the weight of each disjunct
    @variable(s1, lambda[i in lz1,j in lz2]>=0)

    # x=x̄
    @constraint(s1, tu, qu == qubar)
    @constraint(s1, tv, qv== qvbar)

    # Ax+g(y) ⩽0, g2(y)⩽0    
    @constraint(s1, eq1[i in lz1,j in lz2], dot_up[i,j] <= dot_qu[i,j] )
    @constraint(s1, eq2[i in lz1,j in lz2], dot_vp[i,j] <= dot_qv[i,j])
    @NLconstraint(s1, eq3[i in lz1,j in lz2], -( (1-ϵ)* lambda[i,j]  + ϵ)*log(1+dot_u[i,j]/((1-ϵ) *lambda[i,j] + ϵ))+dot_up[i,j] <= 0)
    @NLconstraint(s1, eq4[i in lz1,j in lz2], -( (1-ϵ)* lambda[i,j]  + ϵ)*log(1+dot_v[i,j]/((1-ϵ) *lambda[i,j] + ϵ))+dot_vp[i,j] <= 0)
    @constraint(s1, eq5[i in lz1,j in lz2], dot_u[i,j]<= 10*dot_z1[i,j])
    @constraint(s1, eq6[i in lz1,j in lz2], dot_v[i,j]<= 12*dot_z2[i,j])
    @constraint(s1, eq7[i in lz1,j in lz2], dot_u[i,j]+dot_v[i,j]<=ds*lambda[i,j])

    #x=∑ẋ, y=∑ẏ, ∑λ=1
    @constraint(s1, sum(dot_u[i,j] for i in lz1 for j in lz2) == u)
    @constraint(s1, sum(dot_v[i,j] for i in lz1 for j in lz2)== v)
    @constraint(s1, sum(dot_up[i,j] for i in lz1 for j in lz2)== up)
    @constraint(s1, sum(dot_vp[i,j] for i in lz1 for j in lz2)== vp)
    @constraint(s1, sum(dot_qu[i,j] for i in lz1 for j in lz2)== qu)
    @constraint(s1, sum(dot_qv[i,j] for i in lz1 for j in lz2)==qv)
    @constraint(s1, sum(dot_z1[i,j] for i in lz1 for j in lz2) == z1)
    @constraint(s1, sum(dot_z2[i,j] for i in lz1 for j in lz2 ) == z2 )
    @constraint(s1, sum(lambda[i,j] for i in lz1 for j in lz2)== 1)

    #0⩽ẋ⩽xᵘᵇλ    0⩽ẏ⩽yᵘᵇλ
    @constraint(s1, ub1[i in lz1,j in lz2], dot_up[i,j] <= lambda[i,j] * 3 )
    @constraint(s1, ub2[i in lz1,j in lz2], dot_vp[i,j] <= lambda[i,j] * 3)
    @constraint(s1, ub3[i in lz1,j in lz2], dot_v[i,j] <= lambda[i,j] * 12)
    @constraint(s1, ub4[i in lz1,j in lz2], dot_u[i,j] <= lambda[i,j] * 10)
    @constraint(s1, ub5[i in lz1,j in lz2], dot_qu[i,j] <=lambda[i,j]*2 )
    @constraint(s1, ub6[i in lz1,j in lz2], dot_qv[i,j] <=lambda[i,j]*3)
    @constraint(s1, ub7[i in lz1,j in lz2], dot_z1[i,j] <= lambda[i,j])
    @constraint(s1, ub8[i in lz1,j in lz2], dot_z2[i,j] <= lambda[i,j])

    #yⱼ=0 or 1
    @constraint(s1, spec1[i in lz1, j in lz2], dot_z1[i,j] == (i -1) * lambda[i,j]  )
    @constraint(s1, spec2[i in lz1, j in lz2], dot_z2[i,j] == (j -1) * lambda[i,j]  )

    #objective
    @objective(s1, Min, probs*( 15 * z1 + 12 * z2  + 3* ( u + v) - 50 * (up + vp )))

    return s1 
end 
