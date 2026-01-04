"""
    PtoC(α, κ::AbstractVector{<:Integer}, m::Integer) -> Any

Convert P normalization to C normalization for Jack polynomials.

This computes the ratio between P-normalized and C-normalized Jack polynomials.
Unlike JtoC, this only uses Uhook (not Lhook).

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `m`: Number of variables (can be 0 for symbolic case)

# Examples
```jldoctest
julia> using SymPy
julia> @syms a
julia> MOPS.PtoC(a, [2], 0)
# Returns symbolic expression
```
"""
function PtoC(α, κ::AbstractVector{<:Integer}, m::Integer)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    r = 1
    for i in 1:length(κ)
        for j in 1:κ[i]
            r = r * Uhook(α, κ, j, i)  # Only Uhook, not Lhook
        end
    end
    
    j = sum(κ)
    return r / (α^j * factorial(j))
end

