using JuMP

include("master.jl")
include("ubsub.jl")
include("sub.jl")

m = generate_master()


all_ds = [1.5 2]
sub_problem=[]
ub_problem =[]
ub = 1e3
lb = -1e3
ub_record = []
lb_record = []
mult_p1 = []
mult_p2 = []
mult_x1 = []
mult_x2 = []

g = []
sub_stat=[]
ub_stat=[]

#generate subproblem 
sub_problem = []
ub_problem = []
for i = 1:2
    push!(sub_problem, generate_sub())
end

iter = 1
while ub >= lb * 0.999
    #solve master problem 
    solve(m)
    global p1 = getvalue(getindex(m, :p1))
    global p2 = getvalue(getindex(m, :p2))
    global x1 = getvalue(getindex(m, :x1))
    global x2 = getvalue(getindex(m, :x2))

    #change the first stage decisions and solve subproblem 
    temp_mult_p1 = [0.0,0.0]
    temp_mult_p2 = [0.0,0.0]
    temp_mult_x1 = [0.0, 0.0]
    temp_mult_x2 = [0.0, 0.0]
    temp_mult_g = [0.0,0.0]
    temp_sub_stat = [:ifOptimal, :ifOptimal]
    temp_ub_stat = [:ifOptimal, :ifOptimal]

    for i = 1:2
        sub_problem[i] = generate_sub(p1bar=p1, p2bar=p2, x1bar=x1, x2bar=x2, demand=all_ds[i])
        temp_stat1 =solve(sub_problem[i])
        temp_sub_stat[i] =  temp_stat1
        temp_mult_p1[i] = getdual(getindex(sub_problem[i], :t1))
        temp_mult_p2[i] = getdual(getindex(sub_problem[i], :t2))
        temp_mult_x1[i] = getdual(getindex(sub_problem[i], :t3))
        temp_mult_x2[i] = getdual(getindex(sub_problem[i], :t4))        
        temp_mult_g[i] = getobjectivevalue(sub_problem[i]) - p1 * temp_mult_p1[i] - p2 * temp_mult_p2[i]- x1 * temp_mult_x1[i] - x2 * temp_mult_x2[i]
    end


    if iter == 1
        mult_p1 = temp_mult_p1
        mult_p2 = temp_mult_p2
        mult_x1 = temp_mult_x1
        mult_x2 = temp_mult_x2
        g = temp_mult_g
        sub_stat = temp_sub_stat
        
    else
        mult_p1 = [mult_p1 temp_mult_p1]
        mult_p2 = [mult_p2 temp_mult_p2]
        mult_x1 = [mult_x1 temp_mult_x1]
        mult_x2 = [mult_x2 temp_mult_x2]        
        g = [g temp_mult_g]
        sub_stat = [sub_stat temp_sub_stat]
    end 


    #record objective
    push!(ub_record, getobjectivevalue(sub_problem[1])+getobjectivevalue(sub_problem[2]))
    push!(lb_record, getobjectivevalue(m))
    if lb < lb_record[iter]
        lb = lb_record[iter]
    end
    if ub > ub_record[iter]
        ub = ub_record[iter]
    end
    println("the lower bound  ", lb, "\nthe upper bond is ", ub)

    #generate new master problem
    m = generate_master(mult_p1=mult_p1, mult_p2=mult_p2,mult_x1=mult_x1, mult_x2=mult_x2, g=g, iter=1:iter)
    iter = iter + 1

    # if iter >10
    #     break
    # end
end 
solve(m)
p1 = getvalue(getindex(m, :p1))
p2 = getvalue(getindex(m, :p2))
x1 = getvalue(getindex(m, :x1))
x2 = getvalue(getindex(m, :x2))

for i = 1:2
    push!(ub_problem, generate_ubsub(p1bar=p1, p2bar=p2, x1bar=x1, x2bar=x2, demand=all_ds[i]))
    solve(ub_problem[i])
end
global_ub = getobjectivevalue(ub_problem[1]) + getobjectivevalue(ub_problem[2])

    #solve master problem 

println(sub_stat)
println(lb_record)
println(ub_record)
println("global_ub")
println(global_ub)
println("mult_p1_s1", mult_p1[1,:])
println("mult_p1_s2", mult_p1[2,:])
println("mult_p2_s1", mult_p2[1,:])
println("mult_p2_s2", mult_p2[2,:])
println("mult_x1_s1", mult_x1[1,:])
println("mult_x1_s2", mult_x1[2,:])
println("mult_x2_s1", mult_x2[1,:])
println("mult_x2_s2", mult_x2[2,:])
println("g_s1", g[1,:])
println("g_s2", g[2,:])
println("p1\n", p1)
println("p2\n", p2)
println("x1\n", x1)
println("x2\n", x2)

include("node1.jl")
ub = 1e3 
while ub >= lb * 0.999
    #solve master problem 
    solve(m)
    global p1 = getvalue(getindex(m, :p1))
    global p2 = getvalue(getindex(m, :p2))
    global x1 = getvalue(getindex(m, :x1))
    global x2 = getvalue(getindex(m, :x2))

    #change the first stage decisions and solve subproblem 
    temp_mult_p1 = [0.0,0.0]
    temp_mult_p2 = [0.0,0.0]
    temp_mult_x1 = [0.0, 0.0]
    temp_mult_x2 = [0.0, 0.0]
    temp_mult_g = [0.0,0.0]
    temp_sub_stat = [:ifOptimal, :ifOptimal]
    temp_ub_stat = [:ifOptimal, :ifOptimal]
    println(p1)
    println(p2)
    for i = 1:2
        sub_problem[i] = generate_node1(p1bar=p1, p2bar=p2, x1bar=x1, x2bar=x2, demand=all_ds[i])
        temp_stat1 =solve(sub_problem[i])
        temp_sub_stat[i] =  temp_stat1
        temp_mult_p1[i] = getdual(getindex(sub_problem[i], :t1))
        temp_mult_p2[i] = getdual(getindex(sub_problem[i], :t2))
        temp_mult_x1[i] = getdual(getindex(sub_problem[i], :t3))
        temp_mult_x2[i] = getdual(getindex(sub_problem[i], :t4))        
        temp_mult_g[i] = getobjectivevalue(sub_problem[i]) - p1 * temp_mult_p1[i] - p2 * temp_mult_p2[i]- x1 * temp_mult_x1[i] - x2 * temp_mult_x2[i]
    end



    mult_p1 = [mult_p1 temp_mult_p1]
    mult_p2 = [mult_p2 temp_mult_p2]
    mult_x1 = [mult_x1 temp_mult_x1]
    mult_x2 = [mult_x2 temp_mult_x2]        
    g = [g temp_mult_g]
    sub_stat = [sub_stat temp_sub_stat]



    #record objective
    push!(ub_record, getobjectivevalue(sub_problem[1])+getobjectivevalue(sub_problem[2]))
    push!(lb_record, getobjectivevalue(m))
    if lb < lb_record[iter]
        lb = lb_record[iter]
    end
    if ub > ub_record[iter]
        ub = ub_record[iter]
    end
    println("the lower bound  ", lb, "\nthe upper bond is ", ub)

    #generate new master problem
    m = generate_master(mult_p1=mult_p1, mult_p2=mult_p2,mult_x1=mult_x1, mult_x2=mult_x2, g=g, iter=1:iter,node1=1.0346717044982208)
    iter = iter + 1

    # if iter >10
    #     break
    # end
end 
solve(m)
p1 = getvalue(getindex(m, :p1))
p2 = getvalue(getindex(m, :p2))
x1 = getvalue(getindex(m, :x1))
x2 = getvalue(getindex(m, :x2))

for i = 1:2
    ub_problem[i] = generate_ubsub(p1bar=p1, p2bar=p2, x1bar=x1, x2bar=x2, demand=all_ds[i])
    solve(ub_problem[i])    
end
global_ub = getobjectivevalue(ub_problem[1]) + getobjectivevalue(ub_problem[2])

println(sub_stat)
println(lb_record)
println(ub_record)
println("global_ub")
println(global_ub)
println("mult_p1_s1", mult_p1[1,:])
println("mult_p1_s2", mult_p1[2,:])
println("mult_p2_s1", mult_p2[1,:])
println("mult_p2_s2", mult_p2[2,:])
println("mult_x1_s1", mult_x1[1,:])
println("mult_x1_s2", mult_x1[2,:])
println("mult_x2_s1", mult_x2[1,:])
println("mult_x2_s2", mult_x2[2,:])
println("g_s1", g[1,:])
println("g_s2", g[2,:])
println("p1\n", p1)
println("p2\n", p2)
println("x1\n", x1)
println("x2\n", x2)
