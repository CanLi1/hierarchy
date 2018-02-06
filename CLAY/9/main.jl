

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

push!(node_list, root)
global_ub = root.V_node
global_lb = root.v_node
optimal_x = root.x_opt
optimal_y = root.y_opt


if abs(root.v_node - root.V_node) / abs(root.v_node) < 0.001
    push!(dump_node_list, root)
    node_list = []
end

while length(node_list) >0
    #pick the node with least lower bound
    node_to_branch = 1
    v = 1e8
    for i in 1:length(node_list)    
        if node_list[i].v_node < v 
            v = node_list[i].v_node
            node_to_branch = i 
        end
    end

    newnode1 = deepcopy(node_list[node_to_branch])
    newnode2 = deepcopy(node_list[node_to_branch])
    #select the most fractional variable to branch 
    x_frac = zeros(length(rectangles))
    y_frac = zeros(length(rectangles))
    for i in 1:length(rectangles)
    	x_frac[i] = min(node_list[node_to_branch].xub[i] - node_list[node_to_branch].x_opt[i], node_list[node_to_branch].x_opt[i] - node_list[node_to_branch].xlb[i])
    	y_frac[i] = min(node_list[node_to_branch].yub[i] - node_list[node_to_branch].y_opt[i], node_list[node_to_branch].y_opt[i] - node_list[node_to_branch].ylb[i])
    end

    d_x = get_frac_var(x_frac)
    d_y = get_frac_var(y_frac)

    if d_x[:frac] > d_y[:frac]
        newnode1.xlb[d_x[:index]] = node_list[node_to_branch].x_opt[d_x[:index]]
        newnode2.xub[d_x[:index]] = node_list[node_to_branch].x_opt[d_x[:index]]
    else
        newnode1.ylb[d_y[:index]] = node_list[node_to_branch].y_opt[d_y[:index]]
        newnode2.yub[d_y[:index]] = node_list[node_to_branch].y_opt[d_y[:index]]
    end
    solve_node(newnode1)
    solve_node(newnode2)
    push!(node_list, newnode1)
    push!(node_list, newnode2)
    #update ub 
    if global_ub > newnode1.V_node
        global_ub = newnode1.V_node
        optimal_x = newnode1.x_opt
		optimal_y = newnode1.y_opt
    end
    if global_ub > newnode2.V_node
        global_ub = newnode2.V_node
        optimal_x = newnode2.x_opt
		optimal_y = newnode2.y_opt       
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
global_lb = 1e8 
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
println(optimal_x)
println(optimal_y)
println("\n\n")
for i in 1:length(dump_node_list)
    print("node ")
    println(i)
    println(dump_node_list[i].v_node)
    println(dump_node_list[i].V_node)
    println("number of Benders iterations")
    println(length(dump_node_list[i].mult_x))
    println("variable upper bound ")
    println(dump_node_list[i].xub)
    println(dump_node_list[i].yub)
    println("variable lower bound ")
    println(dump_node_list[i].xub)
    println(dump_node_list[i].xlb)
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






