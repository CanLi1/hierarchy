

addprocs(3)
@everywhere include("input1.jl")
@everywhere include("sub1.jl")
@everywhere include("master.jl")
@everywhere include("ubsub.jl")
@everywhere include("nlprelax.jl")
@everywhere include("util.jl")
@everywhere include("mosek_nlprelax.jl")

#generate subproblem
sub_problem = []
ub_problem = []
sub_stat = []
ub_stat = []
UB = 1e8
LB = -1e8
relax_UB = 1e8
relax_LB = -1e8
djc_scenarios = []
mult_yf = []
mult_n = []
mult_v = []
g = []

sub_obj_record = []
a = now()
b = now()
sub_time = a - b
master_time = a - b
ubsub_time = a - b 
resolve_sub_time = a - b 

for s in scenarios
    push!(djc_scenarios, [[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16],[17],[18],[19],[20],[21],[22],[23],[24]])
end

#generate subproblem and master problem 
for s in scenarios
    push!(sub_problem, generate_nlprelax(Q=Q[:,s], prob =prob[s]))
end

m = generate_master()
obj_master = []
a = now()
while relax_UB > relax_LB + 1
    c = now()
    solve(m)
    d = now()
    master_time = master_time + d - c 
    relax_LB = getobjectivevalue(m)
    push!(obj_master, getobjectivevalue(m))
    yfbarr = getvalue(getindex(m, :yf))
    nbarr = zeros(length(stages))
    for int in integer
        for j in stages
            if yfbarr[int,j] >0.99
                nbarr[j] = log(int)
            end
        end
    end
    vbarr = getvalue(getindex(m, :v))

    temp_mult_yf = []
    temp_mult_n = []
    temp_mult_v = []
    temp_g = []
    temp_sub_stat = []

    #initialize 
    for s in scenarios
    	push!(temp_mult_yf, zeros(length(integer), length(stages)))
    	push!(temp_mult_n, zeros(length(stages)))
        push!(temp_mult_v, zeros(length(stages)))
    	push!(temp_g, 0.0)

    	#update first stage decisions
        for j in stages
            setvalue(getindex(sub_problem[s], :nbar)[j], nbarr[j])
            setvalue(getindex(sub_problem[s], :vbar)[j], vbarr[j])
            for int in integer
                setvalue(getindex(sub_problem[s], :yfbar)[int,j], yfbarr[int,j])
            end
        end
    end
    println(yfbarr)
    println(nbarr)
    println(vbarr)
    c = now()
    results = pmap(psolve, sub_problem)
    d = now()
    sub_time = sub_time + d - c

    for s in scenarios
        sub_problem[s] = results[s][:model]
        if results[s][:status] != :Optimal
            # error("NLP solver converges to an infeasible solution")
            sub_problem[s] = generate_mosek_nlprelax(Q=Q[:,s], prob =prob[s])
            for j in stages
                setvalue(getindex(sub_problem[s], :nbar)[j], nbarr[j])
                setvalue(getindex(sub_problem[s], :vbar)[j], vbarr[j])
                for int in integer
                    setvalue(getindex(sub_problem[s], :yfbar)[int,j], yfbarr[int,j])
                end
            end   
            c = now()
            temp = solve(sub_problem[s])
            d = now()
            sub_time = sub_time + d - c
            results[s][:status] = temp
            results[s][:yf_dual] = getdual(getindex(sub_problem[s], :t1))
            results[s][:v_dual] = getdual(getindex(sub_problem[s], :t2))
            results[s][:n_dual] = getdual(getindex(sub_problem[s], :t3))
            results[s][:objective] = getobjectivevalue(sub_problem[s])

            #revert back to ipopt
            sub_problem[s] = generate_nlprelax(Q=Q[:,s], prob =prob[s])
        end  

        push!(temp_sub_stat, results[s][:status])	
        for j in stages
            temp_mult_n[s][j] = results[s][:n_dual][j]
            temp_mult_v[s][j] = results[s][:v_dual][j]
        	for int in integer
        		temp_mult_yf[s][int, j] = results[s][:yf_dual][int, j]
        	end
        end
        temp_g[s] = results[s][:objective]- sum(sum(temp_mult_yf[s][int,j] * yfbarr[int, j] for int in integer)+ temp_mult_n[s][j] * nbarr[j]+ temp_mult_v[s][j] * vbarr[j] for j in stages)
    end

    #update upper bound 
    if relax_UB > sum(results[s][:objective] for s in scenarios)
        relax_UB = sum(results[s][:objective] for s in scenarios)
    end

    push!(sub_obj_record, sum(results[s][:objective] for s in scenarios))

    #update multipliers
    push!(mult_yf, temp_mult_yf)
    push!(mult_n, temp_mult_n)
    push!(mult_v, temp_mult_v)
    push!(g, temp_g)
    push!(sub_stat, temp_sub_stat)

    m = generate_master(mult_yf=mult_yf, mult_n=mult_n, mult_v=mult_v,g=g, iter=1:length(mult_yf))
    println(sub_obj_record)
    println(obj_master)
    # if length(mult_yf)> 50 && obj_master[length(mult_yf)] - obj_master[length(mult_yf)-1] <1e-3
    # 	break
    # end
    # if length(mult_yf) > 1
    #     break
    # end
    
end

temp_length = length(mult_yf)
println(obj_master)
println(sub_obj_record)


yfbar_for_ub = zeros(length(integer), length(stages))
nbar_for_ub = zeros(length(stages))
vbar_for_ub = zeros(length(stages))

while UB > LB * 1.001
    relax_UB = 1e8
    #generate subproblem and master problem 
    sub_problem=[]
    for s in scenarios
        push!(sub_problem, generate_sub(djc=djc_scenarios[s], Q=Q[:,s], prob =prob[s]))
    end
    m = generate_master(mult_yf=mult_yf, mult_n=mult_n, mult_v=mult_v,g=g, iter=1:length(mult_yf))

    while relax_UB > relax_LB + 1e-1
        c = now()
        solve(m)
        d = now()
        master_time = master_time + d - c 
        relax_LB = getobjectivevalue(m)
        push!(obj_master, getobjectivevalue(m))
        yfbarr = getvalue(getindex(m, :yf))
        nbarr = zeros(length(stages))
        for int in integer
            for j in stages
                if yfbarr[int,j] >0.99
                    nbarr[j] = log(int)
                end
            end
        end
        vbarr = getvalue(getindex(m, :v))

        temp_mult_yf = []
        temp_mult_n = []
        temp_mult_v = []
        temp_g = []
        temp_sub_stat = []

        #initialize 
        for s in scenarios
            push!(temp_mult_yf, zeros(length(integer), length(stages)))
            push!(temp_mult_n, zeros(length(stages)))
            push!(temp_mult_v, zeros(length(stages)))
            push!(temp_g, 0.0)

            #update first stage decisions
            for j in stages
                setvalue(getindex(sub_problem[s], :nbar)[j], nbarr[j])
                setvalue(getindex(sub_problem[s], :vbar)[j], vbarr[j])
                for int in integer
                    setvalue(getindex(sub_problem[s], :yfbar)[int,j], yfbarr[int,j])
                end
            end
        end

        c = now()
        results = pmap(psolve, sub_problem)
        d = now()
        sub_time = sub_time + d - c

        for s in scenarios
            sub_problem[s] = results[s][:model]
            if results[s][:status] != :Optimal
                error("NLP solver converges to an infeasible solution")
            end  

            push!(temp_sub_stat, results[s][:status])   
            for j in stages
                temp_mult_n[s][j] = results[s][:n_dual][j]
                temp_mult_v[s][j] = results[s][:v_dual][j]
                for int in integer
                    temp_mult_yf[s][int, j] = results[s][:yf_dual][int, j]
                end
            end
            temp_g[s] = results[s][:objective]- sum(sum(temp_mult_yf[s][int,j] * yfbarr[int, j] for int in integer)+ temp_mult_n[s][j] * nbarr[j]+ temp_mult_v[s][j] * vbarr[j] for j in stages)
        end

        #update upper bound 
        if relax_UB > sum(results[s][:objective] for s in scenarios)
            relax_UB = sum(results[s][:objective] for s in scenarios)
        end

        push!(sub_obj_record, sum(results[s][:objective] for s in scenarios))

        #update multipliers
        push!(mult_yf, temp_mult_yf)
        push!(mult_n, temp_mult_n)
        push!(mult_v, temp_mult_v)
        push!(g, temp_g)
        push!(sub_stat, temp_sub_stat)

        m = generate_master(mult_yf=mult_yf, mult_n=mult_n, mult_v=mult_v,g=g, iter=1:length(mult_yf))
       
        if relax_UB > relax_LB + 1e-1
            yfbar_for_ub = yfbarr
            nbar_for_ub = nbarr
            vbar_for_ub = vbarr
        end
        if length(mult_yf) > temp_length + 200
            break
        end

        
    end 

    #update lower bound
    LB = relax_LB

    #regenerate the sub_problem and solve again
    sub_problem = []
    for s in scenarios
        push!(sub_problem, generate_sub(djc=djc_scenarios[s], Q=Q[:,s], prob =prob[s]))
        for j in stages
            setvalue(getindex(sub_problem[s], :nbar)[j], nbar_for_ub[j])
            setvalue(getindex(sub_problem[s], :vbar)[j], vbar_for_ub[j])
            for int in integer
                setvalue(getindex(sub_problem[s], :yfbar)[int,j], yfbar_for_ub[int,j])
            end
        end
        c = now()
        temp = solve(sub_problem[s])
        d = now()
        resolve_sub_time = resolve_sub_time + d - c 

    end    
    #now check if the solution is in convex hull 
    all_scenarios_in_convex_hull = true

    #break if the subproblem size is too large
    should_break = false

    for s in scenarios
        λ = getvalue(getindex(sub_problem[s], :dot_λ))
        dot_y = getvalue(getindex(sub_problem[s], :dot_y))
        cur_djc = djc_scenarios[s]
        index_d = 1
        y_over_λ = []
        is_in_convex_hull = false
        for d in cur_djc
            y_over_λ_d = []
            is_djc_d_in_convex_hull = true
            for k in 1:2^length(d)
                y_over_λ_d_k = -1 * ones(length(J1))
                if λ[index_d, k] > 1e-5
                    for i in J1
                        y_over_λ_d_k[i] = dot_y[i, index_d, k] / λ[index_d, k]
                    end
                end
                if  !check_integer(y_over_λ_d_k)
                    is_djc_d_in_convex_hull = false
                end
                push!(y_over_λ_d, y_over_λ_d_k)
            end
            if is_djc_d_in_convex_hull 
                is_in_convex_hull = true
            end
            push!(y_over_λ, y_over_λ_d)
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
            for y_over_λ_d in y_over_λ
                frac_var_index = get_disjunction_frac(y_over_λ_d)
                if frac_var_index[:frac] < frac
                    frac = frac_var_index[:frac]
                    least_frac_index = index_d
                end
                index_d = index_d + 1
            end

            #get least_frac_index and index_disjunction 
            binary_index_j = get_disjunction_frac(y_over_λ[least_frac_index])[:index]
            index_disjunction = 1
            index_d = 1
            for disjunction in cur_djc
                if binary_index_j in disjunction
                    index_disjunction = index_d
                end
                index_d = index_d + 1
            end

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
    #obtain an upper bound
    #generate upper bound subproblem
    ub_sub_problem = []
    for s in scenarios
        push!(ub_sub_problem, generate_ubsub(Q=Q[:,s], prob =prob[s], nbar=nbar_for_ub, vbar=vbar_for_ub, yfbar=yfbar_for_ub))
    end

    #solve upper bound subproblem in parallel
    c = now()
    temp_UB = 0.0
    for s in scenarios
        solve(ub_sub_problem[s])
        temp_UB = temp_UB + getobjectivevalue(ub_sub_problem[s])
    end
    d = now()
    ubsub_time = ubsub_time +  d - c 

    if UB > temp_UB
        UB = temp_UB
    end
    #print 
    println("====================")
    b=now()
    println(b-a)
    println("current upper bound is ")
    println(UB)
    println("current lower bound is ")
    println(LB)
    println(djc_scenarios)
    println(obj_master)
    println(sub_obj_record)
    println("Optimal first stage solution")
    print("nbar=")
    println(nbar_for_ub)
    print("vbar=")
    println(vbar_for_ub)
    println("ubsub_time")
    println(ubsub_time)
    println("sub_time")
    println(sub_time)
    println("master_time")
    println(master_time)
    println("resolve_sub_time")
    println(resolve_sub_time)


    if should_break 
        break
    end
    if length(mult_yf) > temp_length + 300 
        break
    end
   

end
b=now()
println(b-a)
println(temp_length)
println(UB)
println(LB)
println(obj_master)
println(sub_obj_record)
println(djc_scenarios)
println("Optimal first stage solution")
print("nbar=")
println(nbar_for_ub)
print("vbar=")
println(vbar_for_ub)
println("ubsub_time")
println(ubsub_time)
println("sub_time")
println(sub_time)
println("master_time")
println(master_time)
println("resolve_sub_time")
println(resolve_sub_time)

# print("mult_yf=")
# println(mult_yf)
# print("mult_n=")
# println(mult_n)
# println("g")
# println(g)


# println(sub_stat)
# println(mult_yf)
# println(mult_n)
# println(g)
# println(xbar_record)
# println(QEbar_record)














