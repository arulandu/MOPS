"""
    Jacobi(α, κ::AbstractVector{<:Integer}, α1, α2, n, norm) -> Any

Compute multivariate Jacobi polynomials.

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
julia> using SymPy
julia> @syms a a1 a2
julia> MOPS.Jacobi(a, [2], a1, a2, 1)
# Returns symbolic expression
```
"""
function Jacobi(α, κ::AbstractVector{<:Integer}, α1, α2, n, norm=nothing)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition κ: must be non-increasing"))
    end
    
    # Determine normalization convention
    sc = if norm === nothing
        :symC  # Default to symC
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
    
    # Determine number of variables (matching Maple: ks:=n, then if type(n,'list') then ks:=nops(n))
    # Maple allows symbolic n, so we preserve it if it's not a vector
    ks = if n isa AbstractVector
        length(n)
    else
        n  # Can be integer or symbolic
    end
    
    # Get all subpartitions of κ
    sp = subPar(κ)
    println("🔍 Jacobi CALLED | α: $α | κ: $κ | α1: $α1 | α2: $α2 | n: $n | norm: $norm | sc: $sc | ks: $ks")
    println("  → subpartitions: $sp")
    
    result = Sym(0)
    
    # Compute and add coefficients for each subpartition
    for s in sp
        ss = sum(s)
        println("  → Processing subpartition s=$s | ss=$ss")
        
        # Compute the coefficient term:
        # ((-1)^ss * Jacobi_c(α, κ, s, ks, α1, α2) / GSFact(α, α1+(ks-1)/α+1, s) *
        #  Jack(α, s, n, sc) / JackIdentity(α, s, ks)) *
        #  GSFact(α, α1+(ks-1)/α+1, κ) * JackIdentity(α, κ, ks)
        
        jacobi_c_val = Jacobi_c(α, κ, s, ks, α1, α2)
        println("    → jacobi_c_val: $jacobi_c_val")
        
        gsfact_s = GSFact(α, α1 + (ks - 1)/α + 1, s)
        println("    → gsfact_s: $gsfact_s")
        
        jack_s = Jack(α, s, n, sc)
        println("    → jack_s: $jack_s")
        
        jack_identity_s = JackIdentity(α, s, ks)
        println("    → jack_identity_s: $jack_identity_s")
        
        gsfact_k = GSFact(α, α1 + (ks - 1)/α + 1, κ)
        println("    → gsfact_k: $gsfact_k")
        
        jack_identity_k = JackIdentity(α, κ, ks)
        println("    → jack_identity_k: $jack_identity_k")
        
        term = (-1)^ss * (gsfact_k/gsfact_s) * (jack_identity_k/jack_identity_s) * jacobi_c_val * jack_s
        term = simplify(term)
        println("    → term (before simplify): $term")
        println("    → term (after simplify): $(simplify(term))")
        
        result = result + simplify(term)
        println("    → cumulative result: $result")
    end
    
    return result
end

