"""
    evalJack(α, κ::AbstractVector{<:Integer}, x::AbstractVector) -> Any

Evaluate a Jack polynomial at given points.

This function recursively evaluates the Jack polynomial for partition `κ`
at the points in `x`.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `x`: Vector of evaluation points

# Examples
```jldoctest
julia> using SymPy
julia> @syms a x1
julia> MOPS.evalJack(a, [2], [x1])
# Returns symbolic expression
```
"""
function evalJack(α, κ::AbstractVector{<:Integer}, x::AbstractVector)
    if length(κ) > length(x)
        return 0
    end
    
    if isempty(κ)
        return 1
    end
    
    if length(x) == 1
        # Base case: single variable
        # Maple: simplify((x[1]^l[1])*(a^(l[1]-1))*`MOPS/SFact`(1/a+1,l[1]-1))
        return simplify(x[1]^κ[1] * α^(κ[1]-1) * SFact(1/α + 1, κ[1] - 1))
    end
    
    # Recursive case
    # Compute uvalid
    κ_conj = conjugate(κ)
    uvalid = egen(κ_conj)
    
    ls = sum(κ)
    lc = κ_conj
    
    result = 0
    
    for u_conj in uvalid
        u = conjugate(u_conj)
        us = sum(u)
        uc = conjugate(u)
        
        # Compute ab (the coefficient)
        ab = 1
        
        # First loop: over u
        for i in 1:length(u)
            for j in 1:u[i]
                if j <= length(uc) && j <= length(lc) && lc[j] == uc[j]
                    ab = ab / (uc[j] - i + α * (u[i] - j + 1))
                else
                    ab = ab / (uc[j] - i + 1 + α * (u[i] - j))
                end
            end
        end
        
        # Second loop: over l (κ)
        for i in 1:length(κ)
            for j in 1:κ[i]
                if j <= length(uc) && j <= length(lc) && lc[j] == uc[j]
                    ab = ab * (lc[j] - i + α * (κ[i] - j + 1))
                else
                    ab = ab * (lc[j] - i + 1 + α * (κ[i] - j))
                end
            end
        end
        
        # Recursive call
        term = x[end]^(ls - us) * ab * evalJack(α, u, x[1:(end-1)])
        result = result + term
    end
    
    return result
end

