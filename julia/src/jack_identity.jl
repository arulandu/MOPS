"""
    JackIdentity(α, κ::AbstractVector{<:Integer}, m) -> Any

Compute the Jack polynomial normalization constant (identity).

This is the normalization constant used in the Jack polynomial expansion.
Supports both integer and symbolic `m` (matching Maple behavior).

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `m`: Number of variables (can be integer or symbolic)

# Examples
```jldoctest
julia> using SymPy
julia> @syms a
julia> MOPS.JackIdentity(a, [2], 1)
# Returns symbolic expression
```
"""
function JackIdentity(α, κ::AbstractVector{<:Integer}, m)
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
    # Maple: simplify((a^(2*j))*j!*(`MOPS/GSFact`(a,m/a,k))/r)
    # Using SymPy for symbolic computation
    return simplify(α^(2*j) * factorial(j) * GSFact(α, m/α, κ) / r)
end

