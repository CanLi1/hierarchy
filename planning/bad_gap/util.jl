function check_integer(array)
    for item in array
        temp = abs(item % 1)
        if temp > 1e-5 && temp < 1-1e-5
            return false
        end
    end
    return true
end