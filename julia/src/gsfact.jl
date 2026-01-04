"""
    GSFact(α, a, s::AbstractVector{<:Integer})

Compute the generalized shifted factorial over a partition.

The generalized shifted factorial extends the definition of factorials to include partition
parameters. More precisely:
```math
\\operatorname{GSFact}(α, a, s) = \\prod_{i=1}^{\\operatorname{length}(s)} \\operatorname{SFact}\\left(a-\\frac{i-1}{α}, s_i\\right)
```

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `a`: The base value (can be numeric or symbolic)
- `s`: The partition

# Examples
```jldoctest
julia> MOPS.GSFact(1, 3, [3, 2, 1])
360
```
"""
function GSFact(α, a, s::AbstractVector{<:Integer})
    if !parvalid(s)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    result = 1
    for i in 1:length(s)
        result = result * SFact(a - (i - 1) / α, s[i])
    end
    return result
end

