using JuMP

include("full_secondstageCH.jl")
include("fullspace.jl")

# m = generate_fullspace()
# stat = solve(m)
# println("objective value")
# println(getobjectivevalue(m))
# println("solver status")
# println(stat)
# println("First stage varialbes")
# println("q:")
# println(getvalue(getindex(m, :q)))
# println("y:")
# println(getvalue(getindex(m, :y)))
# println("\nSecond stage varialbes")
# println("R:")
# println(getvalue(getindex(m, :R)))
# println("P")
# println(getvalue(getindex(m, :P)))
# println("z")
# println(getvalue(getindex(m, :z)))
# println("δ")
# println(getvalue(getindex(m, :δ)))

list_d = [0.3 0.5  0.6  0.7 0.8 0.9 1 1.1 1.2]
list_ϕ = [120 150 180 210 240 270 300 330 360 390]
list_ψ = [ 1 2 3 4 5 6 7 8 9 10 ]
list_γ = [1 0.9 0.8 1.5 2 3 2.5 3.5 4 0.5]
# list_d = [ 0.8 0.9 1 ]
# list_ϕ = [ 290 300 ]
# list_ψ = [ 8 8.5 9 9.5 10 ]
# list_γ = [1 3]
aγ = [50 50 50 50]
aψ = [1 1 1 1]
aβ =[90 80 100 72]
good_d  = []
good_ϕ = []
good_ψ = []
good_gap = []
good_γ = []
is_z_good_record = []
is_δ_good_record = []
is_q_good_record = []
Qu = [2.3 2.8  2 3.2]
for d in list_d
    for ϕ in list_ϕ
        for ψ in list_ψ
            for γ in list_γ
                m = generate_fullspace(ϕ=ϕ, β=aβ*d, ψ=aψ * ψ, γ= aγ * γ)
                solve(m)
                optvalue = getobjectivevalue(m)
                q = getvalue(getindex(m, :q))
                z = getvalue(getindex(m, :z))
                δ = getvalue(getindex(m, :δ))
                is_q_good = false
                is_z_good = false
                is_δ_good = false
                is_gap_good = false
                for i in 1:3
                    if q[i] > 0.1 && abs(q[i] - Qu[i]) > 0.1
                        is_q_good = true
                    end
                    if z[i,1] - z[i,3] > 0.1
                        is_z_good = true
                    end
                    if δ[i] > 0.1
                        is_δ_good = true
                    end 
                end 
                m2 = generate_fullspaceCH(ϕ=ϕ, β=aβ*d, ψ=aψ * ψ, γ= aγ * γ)
                solve(m2)
                lb = getobjectivevalue(m2)
                if (optvalue - lb) / optvalue > 0.01
                    is_gap_good = true
                end 
                if is_gap_good
                    push!(good_d, d)
                    push!(good_ϕ, ϕ)
                    push!(good_ψ, ψ)
                    push!(good_gap, (optvalue - lb) / optvalue)
                    push!(good_γ, γ)
                    push!(is_q_good_record, is_q_good)
                    push!(is_δ_good_record, is_δ_good)
                    push!(is_z_good_record, is_z_good)
                end 
            end
        end
    end 
end

println(good_d)
println(good_ϕ)
println(good_ψ)
println(good_gap)
println(good_γ)
println(is_q_good_record)
println(is_δ_good_record)
println(is_z_good_record)
