"""
    Uhook(α, κ::AbstractVector{<:Integer}, x::Integer, y::Integer)

Compute the upper hook length of a partition at a square.

The upper hook length is given by: `leg(κ, x, y) + α*(1 + arm(κ, x, y))`

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `x`, `y`: Coordinates of the square

# Examples
```jldoctest
julia> MOPS.Uhook(2, [3, 2, 1], 1, 1)
7
```
"""
function Uhook(α, κ::AbstractVector{<:Integer}, x::Integer, y::Integer)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    return leg(κ, x, y) + α * (1 + arm(κ, x, y))
end

