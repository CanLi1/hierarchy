using JuMP
using Pajarito 
using Ipopt
using CPLEX
using KNITRO
using BARON
function generate_fullspace()
    # s1 =  Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_THREADS=12), cont_solver=IpoptSolver(print_level=0)))
    # s1 = Model(solver=KnitroSolver(KTR_PARAM_PAR_NUMTHREADS=24, KTR_PARAM_MAXTIMECPU=5e4))   
    s1 = Model(solver=BaronSolver(maxtime=5e5, threads=12))
    @variable(s1, p1>=0)
    @variable(s1, p2>=0)
    @variable(s1, x1, Bin)
    @variable(s1, x2, Bin)


    @variable(s1, y1[s in scenarios], Bin)
    @variable(s1, y2[s in scenarios], Bin)
    @variable(s1, q1[s in scenarios]>=0)
    @variable(s1, q2[s in scenarios] >=0)
    @variable(s1, δ[s in scenarios] >= 0)
    #first stage constraints
    @constraint(s1, p1 <= 4 * x1)
    @constraint(s1, p2  <= 2*x2)

    # Ax+g(y) ⩽0, g2(y)⩽0
    @NLconstraint(s1, e1[s in scenarios], (q1[s] -3)^2 + (q2[s] - 2)^2 <= 1 + 16*(1-y1[s]))
    @NLconstraint(s1, e2[s in scenarios], (q1[s]-1)^2 + (q2[s])^2 <=1+16*y1[s] )
    @NLconstraint(s1, e3[s in scenarios], (q1[s])^2 + (q2[s] - 1)^2 <= 1 + 16*(1-y2[s]))
    @NLconstraint(s1, e4[s in scenarios], (q1[s]-4)^2 + (q2[s]-1)^2 <=1+16*y2[s] )
    @constraint(s1, e5[s in scenarios], q1[s]<= p1)
    @constraint(s1, e6[s in scenarios], q2[s]<= p2)
    @constraint(s1, e7[s in scenarios], q1[s] + q2[s] + δ[s]>= d[s])
    @objective(s1, Min, sum(prob[s]* (q1[s]-12*q2[s] + 3*y1[s] -3*y2[s] + 100 * δ[s]) for s in scenarios ) + 3 * x1 +3 * x2 + p1 + p2)
    return s1
end