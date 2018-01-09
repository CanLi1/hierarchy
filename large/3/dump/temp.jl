include("input.jl")
include("sub.jl")
include("ubsub.jl")
include("dump_nlprelax.jl")
total_obj = 0.0
nlprelax_obj = 0.0
cnf_obj = 0.0
obj = []
QE = zeros(2,4)
x = zeros(2,4)
QE[1,1] = 0.5424941419572825
QE[1,2] = 0.16398287866564684
QE[1,3] = 1.6321107729352703e-8
QE[1,4] = 0.03492834168364747
QE[2,1] = 1.6350974519771217e-8
QE[2,2] = 0.7005474899958465
QE[2,3] = 0.6944445176430836
QE[2,4] = 0.0552623726270124
cnf_dual = []
nlprelax_dual = []
for i in 1:2
	for j in 1:4
		if QE[i,j] > 1e-5
			x[i,j] = 1.0
		end
	end
end
QE[2,3] = 0.0944445176430836
QE[2,4] = 0.0052623726270124
model_gen = []
for temp_s in scenarios
	a = now()
	s = generate_ubsub(demand=D[1:2, 1:6, temp_s], price=phi, prob=prob[temp_s], xbar=x, Qbar=QE)
	b = now()
	push!(model_gen, b-a)
	solve(s)
	push!(obj, getobjectivevalue(s))
	total_obj =total_obj + getobjectivevalue(s)

end

for temp_s in scenarios
	a = now()
	s = generate_sub(demand=D[1:2, 1:6, temp_s], price=phi, prob=prob[temp_s], xbar=x, Qbar=QE)
	b = now()
	push!(model_gen, b-a)
	solve(s)
	push!(obj, getobjectivevalue(s))
	cnf_obj =cnf_obj + getobjectivevalue(s)
	push!(cnf_dual,getdual(getindex(s, :t1)))
	push!(cnf_dual,getdual(getindex(s, :t2)))
end

for temp_s in scenarios
	a = now()
	s = generate_nlprelax(demand=D[1:2, 1:6, temp_s], price=phi, prob=prob[temp_s], xbar=x, Qbar=QE)
	b = now()
	push!(model_gen, b-a)
	solve(s)
	push!(obj, getobjectivevalue(s))
	nlprelax_obj =nlprelax_obj + getobjectivevalue(s)
	push!(nlprelax_dual,getdual(getindex(s, :t1)))
	push!(nlprelax_dual,getdual(getindex(s, :t2)))
end
println(obj)
println(total_obj)
println(cnf_obj)
println(nlprelax_obj)
println(cnf_dual)
println(nlprelax_dual)

