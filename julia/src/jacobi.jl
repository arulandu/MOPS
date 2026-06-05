"""
    Jacobi(α, κ::AbstractVector{<:Integer}, α1, α2, n, norm=nothing) -> Any

Compute multivariate Jacobi polynomials (in Jack-C normalization by default).

Computes the Jacobi polynomial for partition `κ` with parameters `α1` and `α2`,
evaluated at points `n` (or with `n` variables if `n` is an integer).

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `α1`: First Jacobi parameter
- `α2`: Second Jacobi parameter
- `n`: Number of variables (integer) or evaluation points (vector)
- `norm`: Normalization convention (optional, defaults to 'symC')
  - Can be 'J', 'P', 'C', 'symJ', 'symP', 'symC', or 'mJ', 'mP', 'mC'

# Examples
```jldoctest
julia> using MOPS
julia> @syms a a1 a2
julia> MOPS.Jacobi(a, [2], a1, a2, 1)
# Returns symbolic expression
```
"""
function Jacobi(α, κ::AbstractVector{<:Integer}, α1, α2, n, norm=nothing)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition κ: must be non-increasing"))
    end
    
    sc = if norm === nothing
        :symC
    else
        norm_str = string(norm)
        if startswith(norm_str, "m")
            # If starts with 'm', extract the rest (e.g., 'mJ' -> 'J')
            Symbol(norm_str[2:end])
        else
            # Otherwise, prepend 'sym' (e.g., 'J' -> 'symJ')
            Symbol("sym" * norm_str)
        end
    end
    
    ks = n isa AbstractVector ? length(n) : n
    sp = subPar(κ)
    result = 0
    gsfact_k = GSFact(α, α1 + (ks - 1)/α + 1, κ)
    jack_identity_k = JackIdentity(α, κ, ks)
    
    for s in sp
        ss = sum(s)
        jacobi_c_val = Jacobi_c(α, κ, s, ks, α1, α2)
        gsfact_s = GSFact(α, α1 + (ks - 1)/α + 1, s)
        jack_s = Jack(α, s, n, sc)
        jack_identity_s = JackIdentity(α, s, ks)
        term = (-1)^ss * (gsfact_k/gsfact_s) * (jack_identity_k/jack_identity_s) * jacobi_c_val * jack_s
        result += term
    end
    
    return expand(result)
end

