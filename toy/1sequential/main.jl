using JuMP

include("master.jl")
include("ubsub.jl")
include("sub.jl")

m = generate_master()


all_ds = [20, 15]
sub_problem=[]
ub_problem =[]
ub = 0.0
lb = 0.0-1e3
ub_record = []
lb_record = []
mult_v = []
mult_u = []
g = []
sub_stat=[]
ub_stat=[]
time_gen = []
#generate subproblem 
sub_problem = []
ub_problem = []
for i = 1:2
    push!(sub_problem, generate_sub(ds=all_ds[i]))
end
    

for i = 1:2
    push!(ub_problem, generate_ubsub(ds=all_ds[i]))
end


iter = 1
while ub >= lb * 0.999
    #solve master problem 
    solve(m)
    qu = getvalue(getindex(m, :qu))
    qv = getvalue(getindex(m, :qv))

    #change the first stage decisions and solve subproblem 
    temp_mult_v = [0.0,0.0]
    temp_mult_u = [0.0,0.0]
    temp_mult_g = [0.0,0.0]
    temp_sub_stat = [:ifOptimal, :ifOptimal]
    temp_ub_stat = [:ifOptimal, :ifOptimal]

    for i = 1:2
        # setvalue(getindex(sub_problem[i], :qub), qu)
        # setvalue(getindex(sub_problem[i], :qvb), qv)
        # setvalue(getindex(ub_problem[i], :qub), qu)
        # setvalue(getindex(ub_problem[i], :qvb), qv)
        a = now()
        sub_problem[i] = generate_sub(qubar=qu, qvbar=qv, ds=all_ds[i])
        ub_problem[i] = generate_ubsub(qubar=qu, qvbar=qv, ds=all_ds[i])
        b= now()
        push!(time_gen, b-a)
        temp_stat1 =solve(sub_problem[i])
        temp_stat2 = solve(ub_problem[i])
        temp_sub_stat[i] =  temp_stat1
        temp_ub_stat[i] = temp_stat2
        temp_mult_v[i] = getdual(getindex(sub_problem[i], :tv))
        temp_mult_u[i] = getdual(getindex(sub_problem[i], :tu ))
        temp_mult_g[i] = getobjectivevalue(sub_problem[i]) - qu * temp_mult_u[i] - qv * temp_mult_v[i]
    end


    if iter == 1
        mult_v = temp_mult_v
        mult_u = temp_mult_u
        g = temp_mult_g
        sub_stat = temp_sub_stat
        ub_stat = temp_ub_stat
    else
        mult_v = [mult_v temp_mult_v]
        mult_u = [mult_u temp_mult_u]
        g = [g temp_mult_g]
        sub_stat = [sub_stat temp_sub_stat]
        ub_stat = [ub_stat temp_sub_stat]
    end 


    #record objective
    push!(ub_record, getobjectivevalue(ub_problem[1])+getobjectivevalue(ub_problem[2]) + getvalue(getindex(m, :first_stage_cost)))
    push!(lb_record, getobjectivevalue(m))
    if lb < lb_record[iter]
        lb = lb_record[iter]
    end
    if ub > ub_record[iter]
        ub = ub_record[iter]
    end
    println("the lower bound  ", lb, "\nthe upper bond is ", ub)

    #generate new master problem
    m = generate_master(mult_u=mult_u, mult_v=mult_v, g=g, iter=1:iter)
    iter = iter + 1

    if iter >10
        break
    end
end 


println(sub_stat)
println(ub_stat)
println(lb_record)
println(ub_record)
println(time_gen)
