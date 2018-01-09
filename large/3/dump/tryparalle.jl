addprocs(2)

@everywhere include("smallmodel.jl")

model_array = []
for i = 1:3
	push!(model_array, generate_model())
end

for i = 1:3
	solve(model_array[i])
	println(getvalue(getindex(model_array[i], :x)))
end
results = pmap(psolve, model_array)
println(results[3][:x])