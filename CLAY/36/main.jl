

addprocs(23)
@everywhere include("input.jl")
@everywhere include("dnf.jl")
@everywhere include("knitro_dnf.jl")
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
end_time = now()
println("total time ")
println(end_time - start_time)
println("sub_time")
println(root.sub_time)
println("master_time")
println(root.master_time)
println("ubsub_time")
println(root.ubsub_time)
println(root.v_node)
println(root.V_node)
println(root.x_opt)
println(root.y_opt)

# push!(node_list, root)
# global_ub = root.V_node
# global_lb = root.v_node
# optimal_sol = zeros(4)
# optimal_sol[1] = root.p1_opt
# optimal_sol[2] = root.p2_opt
# optimal_sol[3] = root.x1_opt
# optimal_sol[4] = root.x2_opt 

# if abs(root.v_node - root.V_node) / abs(root.v_node) < 0.001
#     push!(dump_node_list, root)
#     node_list = []
# end

# while length(node_list) >0
#     #pick the node with least lower bound
#     node_to_branch = 1
#     v = 1e6
#     for i in 1:length(node_list)    
#         if node_list[i].v_node < v 
#             v = node_list[i].v_node
#             node_to_branch = i 
#         end
#     end

#     newnode1 = deepcopy(node_list[node_to_branch])
#     newnode2 = deepcopy(node_list[node_to_branch])
#     #select the most fractional variable to branch 
#     p1_frac = min(node_list[node_to_branch].vars_ub[1] - node_list[node_to_branch].p1_opt, node_list[node_to_branch].p1_opt - node_list[node_to_branch].vars_lb[1])
#     p2_frac = min(node_list[node_to_branch].vars_ub[2] - node_list[node_to_branch].p2_opt, node_list[node_to_branch].p2_opt - node_list[node_to_branch].vars_lb[2])
#     if p1_frac > p2_frac
#         newnode1.vars_lb[1] = node_list[node_to_branch].p1_opt
#         newnode2.vars_ub[1] = node_list[node_to_branch].p1_opt
#     else
#         newnode1.vars_lb[2] = node_list[node_to_branch].p2_opt
#         newnode2.vars_ub[2] = node_list[node_to_branch].p2_opt
#     end
#     solve_node(newnode1)
#     solve_node(newnode2)
#     push!(node_list, newnode1)
#     push!(node_list, newnode2)
#     #update ub 
#     if global_ub > newnode1.V_node
#         global_ub = newnode1.V_node
#         optimal_sol[1] = newnode1.p1_opt
#         optimal_sol[2] = newnode1.p2_opt
#         optimal_sol[3] = newnode1.x1_opt
#         optimal_sol[4] = newnode1.x2_opt 
#     end
#     if global_ub > newnode2.V_node
#         global_ub = newnode2.V_node
#         optimal_sol[1] = newnode2.p1_opt
#         optimal_sol[2] = newnode2.p2_opt
#         optimal_sol[3] = newnode2.x1_opt
#         optimal_sol[4] = newnode2.x2_opt         
#     end

#     node_to_fathom = []
#     push!(node_to_fathom, node_to_branch)
#     node_list[node_to_branch].fathom_type = :bybranch
#     push!(dump_node_list, node_list[node_to_branch])
#     for i in 1:length(node_list)
#         if i != node_to_branch
#             if global_ub < node_list[i].v_node
#                 push!(node_to_fathom, i)
#                 node_list[i].fathom_type = :bybound
#                 push!(dump_node_list, node_list[i])
#             elseif abs(node_list[i].v_node - node_list[i].V_node) / abs(node_list[i].v_node) < 0.001
#                 node_list[i].fathom_type = :byoptimality
#                 push!(node_to_fathom, i)
#                 push!(dump_node_list, node_list[i])
#             end
#         end
#     end

#     new_node_list = []
#     for i in 1:length(node_list)
#         if !(i in node_to_fathom)
#             push!(new_node_list, node_list[i])
#         end
#     end
#     node_list = new_node_list 


# end

# end_time = now()
#     #update lower bound 
# global_lb = 1e6 
# for n in dump_node_list
#     if n.fathom_type == :byoptimality && n.v_node < global_lb
#         global_lb = n.v_node
#     end
# end
# println("==========================================================")
# print("total number of nodes visited is ")
# println(length(dump_node_list))
# print("global upper bound is ")
# println(global_ub)
# print("global lower bound is ")
# println(global_lb)
# println("optimal first stage solution is ")
# println(optimal_sol)
# println("\n\n")
# for i in 1:length(dump_node_list)
#     print("node ")
#     println(i)
#     println(dump_node_list[i].v_node)
#     println(dump_node_list[i].V_node)
#     println("number of Benders iterations")
#     println(length(dump_node_list[i].mult_x1[1, :]))
#     println("variable upper bound ")
#     println(dump_node_list[i].vars_ub)
#     println("variable lower bound ")
#     println(dump_node_list[i].vars_lb)
#     println("fathom_type")
#     println(dump_node_list[i].fathom_type)
#     println("\n\n")
# end

# #solution time 
# println("total running time ")
# println(end_time - start_time)

# for i in 1:length(dump_node_list)
#     master_time = master_time + dump_node_list[i].master_time
#     sub_time = sub_time + dump_node_list[i].sub_time
#     ubsub_time = ubsub_time + dump_node_list[i].ubsub_time
# end 
# println("master_time")
# println(master_time)
# println("sub_time")
# println(sub_time)
# println("ubsub_time")
# println(ubsub_time)






