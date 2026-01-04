"""
    Lhook(α, κ::AbstractVector{<:Integer}, x::Integer, y::Integer)

Compute the lower hook length of a partition at a square.

The lower hook length is given by: `leg(κ, x, y) + 1 + α*(arm(κ, x, y))`

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `x`, `y`: Coordinates of the square

# Examples
```jldoctest
julia> MOPS.Lhook(2, [3, 2, 1], 1, 1)
8
```
"""
function Lhook(α, κ::AbstractVector{<:Integer}, x::Integer, y::Integer)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    return leg(κ, x, y) + 1 + α * arm(κ, x, y)
end

