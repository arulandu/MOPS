"""
    GBC_cont(α, u::AbstractVector{<:Integer}, n::Integer) -> Any

Continuation function for generalized binomial coefficients.

This is an internal function used by `GBC` to compute continuation values.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `u`: The partition
- `n`: The index to decrement

# Examples
```jldoctest
julia> using SymPy
julia> @syms a
julia> MOPS.GBC_cont(a, [2, 1], 1)
# Returns symbolic expression
```
"""
function GBC_cont(α, u::AbstractVector{<:Integer}, n::Integer)
    if !parvalid(u)
        throw(ArgumentError("Invalid partition u: must be non-increasing"))
    end
    if n < 1 || n > length(u)
        throw(ArgumentError("Index n must be between 1 and length(u)"))
    end
    
    # Create l = u with l[n] decremented
    l = copy(u)
    l[n] = l[n] - 1
    
    # Remove trailing zeros if any
    while !isempty(l) && l[end] == 0
        pop!(l)
    end
    
    ab = Sym(1)  # Use Sym(1) for symbolic multiplication
    
    # Compute conjugates
    lc = conjugate(l)
    uc = conjugate(u)
    
    # First loop: multiply
    for i in 1:length(u)
        for j in 1:u[i]
            if j <= length(uc) && j <= length(lc) && lc[j] == uc[j]
                ab = ab * (uc[j] - i + α * (u[i] - j + 1))
            else
                ab = ab * (uc[j] - i + 1 + α * (u[i] - j))
            end
        end
    end
    
    # Second loop: divide
    for i in 1:length(l)
        for j in 1:l[i]
            if j <= length(uc) && j <= length(lc) && lc[j] == uc[j]
                ab = ab / (lc[j] - i + α * (l[i] - j + 1))
            else
                ab = ab / (lc[j] - i + 1 + α * (l[i] - j))
            end
        end
    end
    
    return simplify(ab)
end

