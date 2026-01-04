"""
    parvalid(k) -> Bool

Validates that a partition is in non-increasing (weakly decreasing) order.

A partition is valid if it is empty or if each element is greater than or equal to the next element.
For symbolic partitions, validation is skipped (matching Maple behavior).

# Examples
```jldoctest
julia> MOPS.parvalid([3, 2, 1])
true

julia> MOPS.parvalid([2, 3, 1])
false

julia> MOPS.parvalid([])
true
```
"""
function parvalid(k::AbstractVector)
    if isempty(k)
        return true
    end
    
    # Only validate ordering for integer partitions (matching Maple behavior)
    # Maple's parvalid checks type(k[i],integer) before comparing
    for i in 1:(length(k) - 1)
        if isa(k[i], Integer) && isa(k[i + 1], Integer)
            if k[i] < k[i + 1]
                @warn "partition $k must be ordered from greatest to least"
                return false
            end
        end
    end
    
    return true
end
