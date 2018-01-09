function check_integer(array)
    for item in array
        temp = abs(item % 1)
        if temp > 1e-5 && temp < 1-1e-5
            return false
        end
    end
    return true
end

function psolve(m::JuMP.Model)
	temp = solve(m)
	d = Dict()
	d[:status] = temp
	d[:objective] = getobjectivevalue(m)
	d[:QE_dual] = getdual(getindex(m, :t1))
	d[:x_dual] = getdual(getindex(m, :t2))
	d[:model] = m 
	return d
end
function psolve_sub(m::JuMP.Model)
	temp = solve(m)
	d = Dict()
	d[:status] = temp
	d[:objective] = getobjectivevalue(m)
	d[:QE_dual] = getdual(getindex(m, :t1))
	d[:x_dual] = getdual(getindex(m, :t2))
	d[:model] = m
	return d
end

function ub_psolve(m::JuMP.Model)
	temp = solve(m)
	d = Dict()
	d[:status] = temp
	d[:objective] = getobjectivevalue(m)
	return d
end