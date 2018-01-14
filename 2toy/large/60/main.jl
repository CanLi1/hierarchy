addprocs(23)
@everywhere include("input.jl")
@everywhere include("sub.jl")
@everywhere include("sub_mosek.jl")
@everywhere include("ubsub.jl")
@everywhere include("master.jl")
@everywhere include("util.jl")

a = now()
b=now()
master_time = a-b
sub_time = a-b 
ubsub_time = a-b 
start_time  = now()

node_list = []
dump_node_list = []



root = node()
solve_node(root)
println(root.v_node)
println(root.V_node)
push!(node_list, root)
global_ub = root.V_node
global_lb = root.v_node
optimal_sol = zeros(4)
optimal_sol[1] = root.p1_opt
optimal_sol[2] = root.p2_opt
optimal_sol[3] = root.x1_opt
optimal_sol[4] = root.x2_opt 

if abs(root.v_node - root.V_node) / abs(root.v_node) < 0.001
    push!(dump_node_list, root)
    node_list = []
end

while length(node_list) >0
    #pick the node with least lower bound
    node_to_branch = 1
    v = 1e6
    for i in 1:length(node_list)    
        if node_list[i].v_node < v 
            v = node_list[i].v_node
            node_to_branch = i 
        end
    end

    newnode1 = deepcopy(node_list[node_to_branch])
    newnode2 = deepcopy(node_list[node_to_branch])
    #select the most fractional variable to branch 
    p1_frac = min(node_list[node_to_branch].vars_ub[1] - node_list[node_to_branch].p1_opt, node_list[node_to_branch].p1_opt - node_list[node_to_branch].vars_lb[1])
    p2_frac = min(node_list[node_to_branch].vars_ub[2] - node_list[node_to_branch].p2_opt, node_list[node_to_branch].p2_opt - node_list[node_to_branch].vars_lb[2])
    if p1_frac > p2_frac
        newnode1.vars_lb[1] = node_list[node_to_branch].p1_opt
        newnode2.vars_ub[1] = node_list[node_to_branch].p1_opt
    else
        newnode1.vars_lb[2] = node_list[node_to_branch].p2_opt
        newnode2.vars_ub[2] = node_list[node_to_branch].p2_opt
    end
    solve_node(newnode1)
    solve_node(newnode2)
    push!(node_list, newnode1)
    push!(node_list, newnode2)
    #update ub 
    if global_ub > newnode1.V_node
        global_ub = newnode1.V_node
        optimal_sol[1] = newnode1.p1_opt
        optimal_sol[2] = newnode1.p2_opt
        optimal_sol[3] = newnode1.x1_opt
        optimal_sol[4] = newnode1.x2_opt 
    end
    if global_ub > newnode2.V_node
        global_ub = newnode2.V_node
        optimal_sol[1] = newnode2.p1_opt
        optimal_sol[2] = newnode2.p2_opt
        optimal_sol[3] = newnode2.x1_opt
        optimal_sol[4] = newnode2.x2_opt         
    end

    node_to_fathom = []
    push!(node_to_fathom, node_to_branch)
    node_list[node_to_branch].fathom_type = :bybranch
    push!(dump_node_list, node_list[node_to_branch])
    for i in 1:length(node_list)
        if i != node_to_branch
            if global_ub < node_list[i].v_node
                push!(node_to_fathom, i)
                node_list[i].fathom_type = :bybound
                push!(dump_node_list, node_list[i])
            elseif abs(node_list[i].v_node - node_list[i].V_node) / abs(node_list[i].v_node) < 0.001
                node_list[i].fathom_type = :byoptimality
                push!(node_to_fathom, i)
                push!(dump_node_list, node_list[i])
            end
        end
    end

    new_node_list = []
    for i in 1:length(node_list)
        if !(i in node_to_fathom)
            push!(new_node_list, node_list[i])
        end
    end
    node_list = new_node_list 


end

end_time = now()
    #update lower bound 
global_lb = 1e6 
for n in dump_node_list
    if n.fathom_type == :byoptimality && n.v_node < global_lb
        global_lb = n.v_node
    end
end
println("==========================================================")
print("total number of nodes visited is ")
println(length(dump_node_list))
print("global upper bound is ")
println(global_ub)
print("global lower bound is ")
println(global_lb)
println("optimal first stage solution is ")
println(optimal_sol)
println("\n\n")
for i in 1:length(dump_node_list)
    print("node ")
    println(i)
    println(dump_node_list[i].v_node)
    println(dump_node_list[i].V_node)
    println("number of Benders iterations")
    println(length(dump_node_list[i].mult_x1[1, :]))
    println("variable upper bound ")
    println(dump_node_list[i].vars_ub)
    println("variable lower bound ")
    println(dump_node_list[i].vars_lb)
    println("fathom_type")
    println(dump_node_list[i].fathom_type)
    println("\n\n")
end

#solution time 
println("total running time ")
println(end_time - start_time)

for i in 1:length(dump_node_list)
    master_time = master_time + dump_node_list[i].master_time
    sub_time = sub_time + dump_node_list[i].sub_time
    ubsub_time = ubsub_time + dump_node_list[i].ubsub_time
end 
println("master_time")
println(master_time)
println("sub_time")
println(sub_time)
println("ubsub_time")
println(ubsub_time)






# println("mult_p2_s1", mult_p2[1,:])
# println("mult_p2_s2", mult_p2[2,:])
# println("mult_x1_s1", mult_x1[1,:])
# println("mult_x1_s2", mult_x1[2,:])
# println("mult_x2_s1", mult_x2[1,:])
# println("mult_x2_s2", mult_x2[2,:])
# println("g_s1", g[1,:])
# println("g_s2", g[2,:])
# println("p1\n", p1)
# println("p2\n", p2)
# println("x1\n", x1)
# println("x2\n", x2)
# dot_p1_over_lambda = zeros(Float64, length(scenarios), 4)
# dot_p2_over_lambda = zeros(Float64, length(scenarios), 4)
# for s in scenarios
#     dot_p1 = getvalue(getindex(sub_problem[s], :dot_p1))
#     dot_p2 = getvalue(getindex(sub_problem[s], :dot_p2)) 
#     lambda = getvalue(getindex(sub_problem[s], :dot_λ))
#     for k in 1:4
#         if lambda[k] < 1e-3
#             dot_p1_over_lambda[s,k] = -1
#             dot_p2_over_lambda[s,k] = -1
#         else
#             dot_p1_over_lambda[s,k] = dot_p1[k] / lambda[k]
#             dot_p2_over_lambda[s,k] = dot_p2[k] / lambda[k]
#         end

#     end
# end
# println("dot_p1_over_lambda\n", dot_p1_over_lambda)
# println("dot_p2_over_lambda\n", dot_p2_over_lambda)

# include("node2.jl")
# ub = 1e3 
# while ub >= lb * 0.999
#     #solve master problem 
#     solve(m)
#     global p1 = getvalue(getindex(m, :p1))
#     global p2 = getvalue(getindex(m, :p2))
#     global x1 = getvalue(getindex(m, :x1))
#     global x2 = getvalue(getindex(m, :x2))

#     #change the first stage decisions and solve subproblem 
#     temp_mult_p1 = [0.0,0.0]
#     temp_mult_p2 = [0.0,0.0]
#     temp_mult_x1 = [0.0, 0.0]
#     temp_mult_x2 = [0.0, 0.0]
#     temp_mult_g = [0.0,0.0]
#     temp_sub_stat = [:ifOptimal, :ifOptimal]
#     temp_ub_stat = [:ifOptimal, :ifOptimal]
#     println(p1)
#     println(p2)
#     for i in scenarios
#         sub_problem[i] = generate_node2(p1bar=p1, p2bar=p2, x1bar=x1, x2bar=x2, demand=d[i])
#         temp_stat1 =solve(sub_problem[i])
#         temp_sub_stat[i] =  temp_stat1
#         temp_mult_p1[i] = getdual(getindex(sub_problem[i], :t1))
#         temp_mult_p2[i] = getdual(getindex(sub_problem[i], :t2))
#         temp_mult_x1[i] = getdual(getindex(sub_problem[i], :t3))
#         temp_mult_x2[i] = getdual(getindex(sub_problem[i], :t4))        
#         temp_mult_g[i] = getobjectivevalue(sub_problem[i]) - p1 * temp_mult_p1[i] - p2 * temp_mult_p2[i]- x1 * temp_mult_x1[i] - x2 * temp_mult_x2[i]
#     end



#     mult_p1 = [mult_p1 temp_mult_p1]
#     mult_p2 = [mult_p2 temp_mult_p2]
#     mult_x1 = [mult_x1 temp_mult_x1]
#     mult_x2 = [mult_x2 temp_mult_x2]        
#     g = [g temp_mult_g]
#     sub_stat = [sub_stat temp_sub_stat]



#     #record objective
#     push!(ub_record, getobjectivevalue(sub_problem[1])+getobjectivevalue(sub_problem[2]))
#     push!(lb_record, getobjectivevalue(m))
#     if lb < lb_record[iter]
#         lb = lb_record[iter]
#     end
#     if ub > ub_record[iter]
#         ub = ub_record[iter]
#     end
#     println("the lower bound  ", lb, "\nthe upper bond is ", ub)

#     #generate new master problem
#     m = generate_master(mult_p1=mult_p1, mult_p2=mult_p2,mult_x1=mult_x1, mult_x2=mult_x2, g=g, iter=1:iter,node2=1.0000013295195567)
#     iter = iter + 1

#     # if iter >10
#     #     break
#     # end
# end 
# solve(m)
# p1 = getvalue(getindex(m, :p1))
# p2 = getvalue(getindex(m, :p2))
# x1 = getvalue(getindex(m, :x1))
# x2 = getvalue(getindex(m, :x2))

# for i in scenarios
#     ub_problem[i] = generate_ubsub(p1bar=p1, p2bar=p2, x1bar=x1, x2bar=x2, demand=d[i])
#     solve(ub_problem[i])    
# end
# global_ub = getobjectivevalue(ub_problem[1]) + getobjectivevalue(ub_problem[2])

# println(sub_stat)
# println(lb_record)
# println(ub_record)
# println("global_ub")
# println(global_ub)
# println("mult_p1_s1", mult_p1[1,:])
# println("mult_p1_s2", mult_p1[2,:])
# println("mult_p2_s1", mult_p2[1,:])
# println("mult_p2_s2", mult_p2[2,:])
# println("mult_x1_s1", mult_x1[1,:])
# println("mult_x1_s2", mult_x1[2,:])
# println("mult_x2_s1", mult_x2[1,:])
# println("mult_x2_s2", mult_x2[2,:])
# println("g_s1", g[1,:])
# println("g_s2", g[2,:])
# println("p1\n", p1)
# println("p2\n", p2)
# println("x1\n", x1)
# println("x2\n", x2)
