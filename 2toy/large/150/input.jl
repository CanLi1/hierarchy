scenarios = 1:150
unit_length = 3.0 / length(scenarios)
prob = ones(length(scenarios)) / length(scenarios)
d = zeros(length(scenarios))
d[1] = unit_length
for s in scenarios
	if s > 1
		d[s] = d[s-1] + unit_length
	end
end
