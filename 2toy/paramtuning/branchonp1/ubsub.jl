using JuMP
using Pajarito
using CPLEX 
using Ipopt

function generate_ubsub(; p1bar=0, p2bar=0, x1bar=0, x2bar=0, prob=0.5, demand=2)
    s1 =  Model(solver=PajaritoSolver(rel_gap=0.00001, mip_solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_THREADS=12), cont_solver=IpoptSolver()))

    @variable(s1, p1>=0)
    @variable(s1, p2>=0)
    @variable(s1, x1)
    @variable(s1, x2)

    @variable(s1, y1, Bin)
    @variable(s1, y2, Bin)
    @variable(s1, q1>=0)
    @variable(s1, q2 >=0)
    @variable(s1, δ>=0)

    @constraint(s1, t1, p1== p1bar)
    @constraint(s1, t2, p2== p2bar)
    @constraint(s1, t3, x1 == x1bar)
    @constraint(s1, t4, x2 == x2bar)

    # Ax+g(y) ⩽0, g2(y)⩽0
    @NLconstraint(s1, e1, (q1 -3)^2 + (q2 - 2)^2 <= 1 + 16*(1-y1))
    @NLconstraint(s1, e2, (q1-1)^2 + (q2 )^2 <=1+16*y1 )
    @NLconstraint(s1, e3, (q1)^2 + (q2 - 1)^2 <= 1 + 16*(1-y2))
    @NLconstraint(s1, e4, (q1-4)^2 + (q2-1)^2 <=1+16*y2 )
    @constraint(s1, e5, q1<= p1)
    @constraint(s1, e6, q2<= p2)
    @constraint(s1, e7, q1 + q2 +δ >= demand)
    @objective(s1, Min, prob* (q1-12*q2 + 3*y1 -3*y2 + 3 * x1 +3 * x2 + p1 + p2+100*δ))
    return s1
end