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
z1_dot_record = []
z2_dot_record = []
qu_dot_record = []
qv_dot_record = []
lamda_record = [] 
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
    temp_dot_z1 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    temp_dot_z2 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    temp_dot_qu = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    temp_dot_qv = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    temp_lambda = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

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
        #get value of zdot 
        z1_dot = getvalue(getindex(sub_problem[i], :dot_z1))
        z2_dot = getvalue(getindex(sub_problem[i], :dot_z2))
        qu_dot = getvalue(getindex(sub_problem[i], :dot_qu))
        qv_dot = getvalue(getindex(sub_problem[i], :dot_qv))        
        ll = getvalue(getindex(sub_problem[i], :lambda))
        for j in 1:2
            for k in 1:2
                index = 4*(i-1) + 2 *(j-1) + k
                temp_dot_z1[index] = z1_dot[k,j]
                temp_dot_z2[index] = z2_dot[k,j]
                temp_dot_qu[index] = qu_dot[k,j]
                temp_dot_qv[index] = qv_dot[k,j]                
                temp_lambda[index] = ll[k,j]
            end
        end
    end


    if iter == 1
        mult_v = temp_mult_v
        mult_u = temp_mult_u
        g = temp_mult_g
        sub_stat = temp_sub_stat
        ub_stat = temp_ub_stat
        z1_dot_record = temp_dot_z1
        z2_dot_record = temp_dot_z2
        qu_dot_record = temp_dot_qu
        qv_dot_record = temp_dot_qv
        lamda_record = temp_lambda
    else
        mult_v = [mult_v temp_mult_v]
        mult_u = [mult_u temp_mult_u]
        g = [g temp_mult_g]
        sub_stat = [sub_stat temp_sub_stat]
        ub_stat = [ub_stat temp_sub_stat]
        z1_dot_record = [z1_dot_record; temp_dot_z1]
        z2_dot_record = [z2_dot_record ;temp_dot_z2]
        qu_dot_record = [qu_dot_record; temp_dot_qu]
        qv_dot_record = [qv_dot_record; temp_dot_qv]
        lamda_record = [lamda_record ;temp_lambda]
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
println("lb record\n", lb_record)
println("ub_record\n",ub_record)
println("z1_dot_record\n",z1_dot_record)
println("z2_dot_record\n",z2_dot_record)
println("lamda_record\n",lamda_record)
println("qu_dot_record\n",qu_dot_record)
println("qv_dot_record\n",qv_dot_record)
println("mult_v\n", mult_v)
println("mult_u\n", mult_u)
println("g\n", g)

