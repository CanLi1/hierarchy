using JuMP
using Mosek

s1 =  Model(solver=MosekSolver(MSK_IPAR_NUM_THREADS=0))
@variable(s1, x>=0)
@variable(s1, y>=0)
@NLconstraint(s1, log(1+x)>=y)
@constraint(s1, y>=2)
@objective(s1, Min, 3*x+6*y)
solve(s1)