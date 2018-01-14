using JuMP

include("master.jl")
include("ubsub.jl")
include("sub.jl")
include("util.jl")
m = generate_master()

prob = [0.3 0.4 0.3]
demand = [10 8.5  6]

#generate subproblem
sub_problem = []
ub_problem = []
sub_stat = []
ub_stat = []
UB = 1e6
LB = -1e6
djc_scenarios = []
mult_q = []
mult_y = []
g = []

for s in 1:3
    push!(djc_scenarios, [[1], [2], [3], [4]])
end


while UB > LB + 1e-3
    while true
        relax_LB = -1e6
        relax_UB = 1e6
        ub_record = []
        lb_record = []
        sub_record = []
        while relax_UB > relax_LB + 0.5 
            obj_master = 0
            solve(m)
            obj_master = getobjectivevalue(m)
            if relax_LB < obj_master
                relax_LB = obj_master
            end
            push!(lb_record, relax_LB)
            global qbarr, ybarr
            qbarr = getvalue(getindex(m, :q))
            ybarr = getvalue(getindex(m, :y))
            temp_mult_q = []
            temp_mult_y = []
            temp_g = []
            global sub_problem=[]
            #initialize
            for s in 1:3
                push!(temp_mult_q, [0.0, 0.0, 0.0, 0.0])
                push!(temp_mult_y, [0.0, 0.0, 0.0, 0.0])
                push!(temp_g, 0)
            end
            temp_sub_stat = []
            obj_sub = 0 
            for s = 1:3
                # println("well before generate sub")
                push!(sub_problem, generate_sub(djc = djc_scenarios[s], demand=demand[s], prob=prob[s], qbar=qbarr, ybar=ybarr))
                temp = solve(sub_problem[s])
                if temp != :Optimal
                    println(ub_record)
                    println(lb_record)
                    println(sub_stat)
                    println(m)
                    println(qbarr)
                    println(ybarr)                    
                    error("Mosek converge to an infeasible solution")
                end
                push!(temp_sub_stat, temp)
                q_dual = getdual(getindex(sub_problem[s], :tq))
                y_dual = getdual(getindex(sub_problem[s], :ty))
                for ii in 1:4
                    temp_mult_q[s][ii] = q_dual[ii]
                    temp_mult_y[s][ii] = y_dual[ii]
                end 
                temp_g[s] = getobjectivevalue(sub_problem[s])-sum(temp_mult_q[s][ii] * qbarr[ii] + temp_mult_y[s][ii] * ybarr[ii] for ii in 1:4)
                obj_sub += getobjectivevalue(sub_problem[s])
            end
            #update multipliers
            push!(mult_q, temp_mult_q)
            push!(mult_y, temp_mult_y)
            push!(g, temp_g)
            push!(sub_stat, temp_sub_stat)

            if obj_sub < relax_UB
                relax_UB = obj_sub
            end
            push!(sub_record, obj_sub)
            push!(ub_record, relax_UB)
            m = generate_master(mult_q=mult_q, mult_y=mult_y, g=g, iter=1:length(mult_q))
            # if length(mult_q) > 2
            #     break
            # end
        end
        println("ub_record", ub_record)
        println("lb_record",lb_record)
        println("sub_stat",sub_stat)
        println("the optimal solution is ")
        println("q ", qbarr)
        println("y ", ybarr)
        all_scenarios_in_convex_hull = true
        for s in 1:3
            println("scenario ", s)
            println("z ", getvalue(getindex(sub_problem[s], :z)))
            # println("dot_z ", getvalue(getindex(sub_problem[s], :dot_z)))
            # println("dot_q ", getvalue(getindex(sub_problem[s], :dot_q)))
            # println("dot_y ", getvalue(getindex(sub_problem[s], :dot_y)))
            # println("dot_λ ", getvalue(getindex(sub_problem[s], :dot_λ)))
            #check if there the z satisfy the integrality constraint
            λ = getvalue(getindex(sub_problem[s], :dot_λ))
            dot_z = getvalue(getindex(sub_problem[s], :dot_z))
            cur_djc = djc_scenarios[s]
            index_d = 1
            z_over_λ = []
            is_in_convex_hull = false
            for d in cur_djc
                z_over_λ_d = []
                is_djc_d_in_convex_hull = true
                for k in 1:2^length(d)
                    z_over_λ_d_k = [-1.0, -1.0, -1.0, -1.0]
                    if λ[index_d, k] > 1e-4
                        for i in 1:4
                            z_over_λ_d_k[i] = dot_z[i, index_d, k] / λ[index_d, k]
                        end
                    end
                    if  !check_integer(z_over_λ_d_k)
                        is_djc_d_in_convex_hull = false
                    end
                    push!(z_over_λ_d, z_over_λ_d_k)
                end
                if is_djc_d_in_convex_hull 
                    is_in_convex_hull = true
                end
                push!(z_over_λ, z_over_λ_d)
                index_d = index_d + 1
            end
            if !is_in_convex_hull
                all_scenarios_in_convex_hull = false
            end
            #if not in convex hull, apply one basic step accroding to heuristic 1
            if !is_in_convex_hull && length(cur_djc) >1
                index_d = 1
                least_frac_index = 1
                frac = 1
                #find the least fractional disjunction
                for z_over_λ_d in z_over_λ
                    frac_var_index = get_disjunction_frac(z_over_λ_d)
                    if frac_var_index[:frac] < frac
                        frac = frac_var_index[:frac]
                        least_frac_index = index_d
                    end
                    index_d = index_d + 1
                end

                #get least_frac_index and index_disjunction 
                binary_index_j = get_disjunction_frac(z_over_λ[least_frac_index])[:index]
                index_disjunction = 1
                index_d = 1
                for disjunction in cur_djc
                    if binary_index_j in disjunction
                        index_disjunction = index_d
                    end
                    index_d = index_d + 1
                end

                println(binary_index_j)
                println(least_frac_index)
                println(index_disjunction)
                println(cur_djc)
                println(z_over_λ)

                #merge least_frac_index and index_disjunction 
                temp_djc = []
                for k in 1:(length(cur_djc)-1)
                    push!(temp_djc, [0])
                end
                temp_djc[1] = [cur_djc[least_frac_index]; cur_djc[index_disjunction]]
                
                #break if size is too large
                # if length(temp_djc[1]) >= 4
                #     should_break = true
                # end
                index_k = 2
                for k in 1:(length(cur_djc))
                    if k != least_frac_index && k != index_disjunction
                        temp_djc[index_k] = cur_djc[k]
                        index_k = index_k + 1
                    end
                end

                djc_scenarios[s] = temp_djc
            end

        end

           
        
        println("*****")
        println(djc_scenarios)
      
        if all_scenarios_in_convex_hull 
            #update  global upper bound and lower bound
            LB = relax_LB
            ub_sub_problem = []
            temp_UB = 0
            for s in 1:3
                push!(ub_sub_problem, generate_ubsub(demand=demand[s], prob=prob[s], qbar=qbarr, ybar=ybarr))
                solve(ub_sub_problem[s])
                temp_UB += getobjectivevalue(ub_sub_problem[s])
            end
            if temp_UB < UB 
                UB = temp_UB
            end
            println("gobal lower bound is ", LB)
            println("global upper bound is ", UB)
            println("the optimal solution is ")
            println("q ", qbarr)
            println("y ", ybarr)
            break
        end
        
    end
    break
end
    

# for s = 1:3
#     push!(ub_problem, generate_ubsub(demand=demand[s], prob=prob[s],qbar=[2.0797 2.227 1.90697 2.38774], ybar=[1 1 1 1]))
#     temp = solve(ub_problem[s])
#     push!(ub_stat, temp)
#     obj_ub += getobjectivevalue(ub_problem[s])
# end


