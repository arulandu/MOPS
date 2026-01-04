"""
    partition_le(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer}) -> Bool

Check if partition `k` is less than or equal to partition `l` in dominance order.

Partition dominance: `k <= l` if for all prefixes, the sum of `k` is <= sum of `l`.

# Examples
```jldoctest
julia> MOPS.partition_le([2, 1], [3, 1])
true

julia> MOPS.partition_le([3, 1], [2, 1])
false
```
"""
function partition_le(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer})
    if !parvalid(k)
        throw(ArgumentError("Invalid partition k: must be non-increasing"))
    end
    if !parvalid(l)
        throw(ArgumentError("Invalid partition l: must be non-increasing"))
    end
    
    min_len = min(length(k), length(l))
    for i in 1:min_len
        if k[i] < l[i]
            return true
        end
        if l[i] < k[i]
            return false
        end
    end
    
    return length(k) <= length(l)
end

"""
    partition_lt(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer}) -> Bool

Check if partition `k` is strictly less than partition `l` in dominance order.
"""
function partition_lt(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer})
    if !parvalid(k)
        throw(ArgumentError("Invalid partition k: must be non-increasing"))
    end
    if !parvalid(l)
        throw(ArgumentError("Invalid partition l: must be non-increasing"))
    end
    
    min_len = min(length(k), length(l))
    for i in 1:min_len
        if k[i] < l[i]
            return true
        end
        if l[i] < k[i]
            return false
        end
    end
    
    return length(k) < length(l)
end

"""
    partition_ge(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer}) -> Bool

Check if partition `k` is greater than or equal to partition `l` in dominance order.
"""
function partition_ge(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer})
    return !partition_lt(k, l)
end

"""
    partition_gt(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer}) -> Bool

Check if partition `k` is strictly greater than partition `l` in dominance order.
"""
function partition_gt(k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer})
    return !partition_le(k, l)
end

