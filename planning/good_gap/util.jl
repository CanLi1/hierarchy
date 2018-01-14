function check_integer(array)
    for item in array
        temp = abs(item % 1)
        if temp > 1e-4 && temp < 1-1e-4
            return false
        end
    end
    return true
end

#return the most fractional variable's fraction in an array and its index 
function get_frac_var(var_array)
	most_frac_index = 1
	frac = 0.0
	index = 1
	for var in var_array
		if var < 0.0
			continue
		end
		if var < 0.50
			if var > frac
				frac = var 
				most_frac_index = index
			end
		else
			if (1-var) > frac
				frac = 1 - var 
				most_frac_index = index
			end
		end
		index = index + 1
	end
	d = Dict()
	d[:frac] = frac
	d[:index] = most_frac_index
	return d
end

function get_disjunction_frac(djc)
	d = Dict()
	d[:frac] = 0
	d[:index] = 1
	for disjunct_y in djc 
		temp = get_frac_var(disjunct_y)
		if temp[:frac] > d[:frac]
			d[:frac] = temp[:frac]
			d[:index] = temp[:index]
		end
	end
	return d
end