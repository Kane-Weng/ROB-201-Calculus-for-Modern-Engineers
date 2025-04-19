function ustep(x)
    value = 1.0
    if x < 0
        value = 0.0
    end
    return value
end

function rect(x)
    value = 0.0
    if abs(x) <= 0.5
        value = 1.0
    end
    return value
end

function tri(x)
    value = 0.0
    if abs(x) < 1
        value = 1.0 - abs(x)
    end
    return value
end