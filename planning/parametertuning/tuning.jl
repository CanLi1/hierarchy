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
list_ϕ = [190 180 170 160 200  210 220 230 240  250 260 270 280 290 300 320 330 340 350 380 400]
list_ψ = [ 8 8.5 9 9.5 10 ]
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
                if is_q_good && is_z_good && is_δ_good
                    push!(good_d, d)
                    push!(good_ϕ, ϕ)
                    push!(good_ψ, ψ)
                    push!(good_gap, (optvalue - lb) / optvalue)
                    push!(good_γ, γ)
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

