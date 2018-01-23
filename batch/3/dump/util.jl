function max(d::Array{Float64,1})
	max_item = d[1]
	for item in d
		if item > max_item
			max_item = item
		end
	end
	return max_item
end

function min(d::Array{Float64,1})
	min_item = d[1]
	for item in d
		if item < min_item
			min_item = item
		end
	end
	return min_item
end