using JuMP
using Pajarito
using CPLEX 
using Ipopt

prob = [0.5 0.5]
qubar = 0
qvbar = 0
ds = 20
function generate_ubsub(;qubar=0, qvbar=0, probs=0.5, ds=20)
    #set up model
    s1 =  Model(solver=PajaritoSolver(mip_solver=CplexSolver(CPX_PARAM_SCRIND=0), cont_solver=IpoptSolver(print_level=0)))

    #continous vars
    @variable(s1, u>=0)
    @variable(s1, v>=0)
    @variable(s1, up>=0)
    @variable(s1, vp>=0)
    @variable(s1, qu)
    @variable(s1, qv)
    
    #binary vars
    @variable(s1, z1, Bin)
    @variable(s1, z2, Bin)

    #transfer the first stage decisions constraints
    @constraint(s1, qu == qubar)
    @constraint(s1, qv== qvbar)

    #second stage constraints
    @constraint(s1, up <= qu )
    @constraint(s1, vp <= qv)
    @NLconstraint(s1, -log(1+u)+up <= 0)
    @NLconstraint(s1, -log(1+v)+vp <=0)
    @constraint(s1, u<= 10*z1)
    @constraint(s1, v<=12*z2)
    @constraint(s1, u+v<=ds)

    #objective
    @objective(s1, Min, probs*( 15 * z1 + 12 * z2  + 3* ( u + v) - 50 * (up + vp )))

    return s1 
end 
