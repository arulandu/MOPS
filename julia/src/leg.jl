"""
    leg(κ, a::Integer, b::Integer)

Compute the leg length of a partition at a square.

The leg length is computed using the conjugate partition and the arm function.

Supports both integer and symbolic partitions (matching Maple behavior).

# Examples
```jldoctest
julia> MOPS.leg([4, 3, 1, 1], 1, 1)
3
```
"""
function leg(κ::AbstractVector, a::Integer, b::Integer)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    κ_conj = conjugate(κ)
    return arm(κ_conj, b, a)
end

