"""
    JtoC(α, κ::AbstractVector{<:Integer}, m::Integer) -> Any

Convert J normalization to C normalization for Jack polynomials.

This computes the ratio between J-normalized and C-normalized Jack polynomials.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `m`: Number of variables (can be 0 for symbolic case)

# Examples
```jldoctest
julia> using SymPy
julia> @syms a
julia> MOPS.JtoC(a, [2], 0)
# Returns symbolic expression
```
"""
function JtoC(α, κ::AbstractVector{<:Integer}, m::Integer)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    r = 1
    for i in 1:length(κ)
        for j in 1:κ[i]
            r = r * Uhook(α, κ, j, i) * Lhook(α, κ, j, i)
        end
    end
    
    j = sum(κ)
    return r / (α^j * factorial(j))
end

