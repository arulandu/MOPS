"""
    arm(κ, a::Integer, b::Integer)

Compute the arm length of a partition at a square.

Given the diagram of the partition `κ`, the arm length of the (a, b) square is the
number of squares to the right of that square, i.e., `arm(κ, a, b) = κ[b] - a`.

Supports both integer and symbolic partitions (matching Maple behavior).

# Examples
```jldoctest
julia> MOPS.arm([4, 3, 1, 1], 2, 1)
2
```
"""
function arm(κ::AbstractVector, a::Integer, b::Integer)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    if b < 1 || b > length(κ)
        throw(BoundsError(κ, b))
    end
    
    return κ[b] - a
end

