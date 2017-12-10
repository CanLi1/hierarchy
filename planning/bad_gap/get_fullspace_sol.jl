using JuMP

include("full_secondstageCH.jl")

m = generate_fullspace()
stat = solve(m)
println("objective value")
println(getobjectivevalue(m))
println("solver status")
println(stat)
println("First stage varialbes")
println("q:")
println(getvalue(getindex(m, :q)))
println("y:")
println(getvalue(getindex(m, :y)))
println("\nSecond stage varialbes")
println("R:")
println(getvalue(getindex(m, :R)))
println("P")
println(getvalue(getindex(m, :P)))
println("z")
println(getvalue(getindex(m, :z)))
println("δ")
println(getvalue(getindex(m, :δ)))

# list_d = [0.1 0.2 0.3 0.5  0.6  0.7 0.8 0.9 1]
# list_ϕ = [200  230  250 270 290 330 350 380 400]
# list_ψ = [1 2 3 4 5 6 7 8 9 10 0.5 0.8 0.3]
# good_d  = []
# good_ϕ = []
# good_ψ = []
# for d in list_d
#     for ϕ in list_ϕ
#         for ψ in list_ψ
#             m = generate_fullspace(ϕ=ϕ, β=aβ*d, ψ=aψ * ψ)
#             solve(m)
#             q = getvalue(getindex(m, :q))
#             z = getvalue(getindex(m, :z))
#             δ = getvalue(getindex(m, :δ))
#             is_q_good = false
#             is_z_good = false
#             is_δ_good = false
#             for i in 1:3
#                 if q[i] > 0.1 && abs(q[i] - Qu[i]) > 0.1
#                     is_q_good = true
#                 end
#                 if z[i,1] - z[i,3] > 0.1
#                     is_z_good = true
#                 end
#                 if δ[i] > 0.1
#                     is_δ_good = true
#                 end 
#             end 
#             if is_q_good && is_z_good && is_δ_good
#                 push!(good_d, d)
#                 push!(good_ϕ, ϕ)
#                 push!(good_ψ, ψ)
#             end 
#         end
#     end 
# end

# println(good_d)
# println(good_ϕ)
# println(good_ψ)

# Any[0.3, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 0.9, 1.0]
# Any[200, 230, 270, 290, 270, 290, 290, 330, 350, 380]
# Any[9.0, 9.0, 8.0, 9.0, 9.0, 8.0, 9.0, 9.0, 9.0, 9.0]