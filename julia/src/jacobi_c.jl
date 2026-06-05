"""
    Jacobi_c(α, κ::AbstractVector{<:Integer}, s::AbstractVector{<:Integer}, m, α1, α2) -> Any

Compute the Jacobi expansion coefficient ``c_{κ,s}`` appearing in

``J_κ = (a)_κ C_κ(I_m) * Σ_{s ⊆ κ} (-1)^|s| c_{κ,s}/(a)_s * C_s(x)/C_s(I_m)``
where ``a = α1 + (m-1)/α + 1``.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The larger partition
- `s`: The smaller partition (subpartition of κ)
- `m`: Number of variables
- `α1`: First Jacobi parameter
- `α2`: Second Jacobi parameter

# Examples
```jldoctest
julia> using MOPS
julia> @syms a a1 a2
julia> MOPS.Jacobi_c(a, [2], [1], 1, a1, a2)
# Returns symbolic expression
```
"""
function Jacobi_c(α, κ::AbstractVector{<:Integer}, s::AbstractVector{<:Integer}, m, α1, α2)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition κ: must be non-increasing"))
    end
    if !parvalid(s)
        throw(ArgumentError("Invalid partition s: must be non-increasing"))
    end
    if !subpar_check(s, κ)
        return 0
    end
    
    if κ == s
        return 1
    end
    
    result = 0
    
    for i in 1:min(length(s) + 1, length(κ))
        if i == length(s) + 1
            si = [s..., 1]
        else
            si = copy(s)
            si[i] = si[i] + 1
        end
        
        # Keep only valid Young diagrams s^(i) contained in κ.
        if si == sort(si, rev=true)
            si_clean = filter(x -> x > 0, si)
            if parvalid(si_clean) && subpar_check(si_clean, κ)
                result += GBC_cont(α, si_clean, i) * Jacobi_c(α, κ, si_clean, m, α1, α2)
            end
        end
    end
    
    ks = sum(κ)
    ss = sum(s)
    val = (α1 + α2 + 2*(m - 1)/α + 2) * (ks - ss) + rho(α, κ) - rho(α, s)
    
    if symbolic_zero(val)
        return 0
    end
    
    return simplify(result / val)
end

