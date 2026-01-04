"""
    isdominate(k::AbstractVector{<:Integer}, s::AbstractVector{<:Integer}) -> Bool

Check if partition `k` dominates partition `s` in dominance order.

A partition `k` dominates `s` if for all prefixes, the sum of `k` is >= the sum of `s`.

# Arguments
- `k`: The first partition
- `s`: The second partition

# Examples
```jldoctest
julia> MOPS.isdominate([3, 1], [2, 2])
true

julia> MOPS.isdominate([2, 2], [3, 1])
false
```
"""
function isdominate(k::AbstractVector{<:Integer}, s::AbstractVector{<:Integer})
    if !parvalid(k)
        throw(ArgumentError("Invalid partition k: must be non-increasing"))
    end
    if !parvalid(s)
        throw(ArgumentError("Invalid partition s: must be non-increasing"))
    end
    
    s1 = 0
    s2 = 0
    max_len = max(length(k), length(s))
    
    for i in 1:max_len
        if i <= length(k)
            s1 += k[i]
        end
        if i <= length(s)
            s2 += s[i]
        end
        if s1 < s2
            return false
        end
    end
    
    return true
end

