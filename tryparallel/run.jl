addprocs(3) #say I have 4 processors available
@everywhere include("ParallelExample.jl")

nmodels = 10

# generate a JuMP.Model array with 50 elements of single-scenario problems
input1 = [5 for i in 1:nmodels]
input2 = [10 for i in 1:nmodels]
input3 = [1 for i in 1:nmodels]

# parallel version
pmap(generate_then_solvemodel, input1, input2, input3)