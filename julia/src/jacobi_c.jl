"""
    Jacobi_c(α, κ::AbstractVector{<:Integer}, s::AbstractVector{<:Integer}, m, α1, α2) -> Any

Compute Jacobi polynomial expansion coefficients.

This is a recursive function that computes the coefficient of the Jack polynomial
for partition `s` in the expansion of the Jacobi polynomial for partition `κ`.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The larger partition
- `s`: The smaller partition (subpartition of κ)
- `m`: Number of variables
- `α1`: First Jacobi parameter
- `α2`: Second Jacobi parameter

# Examples
```jldoctest
julia> using SymPy
julia> @syms a a1 a2
julia> MOPS.Jacobi_c(a, [2], [1], 1, a1, a2)
# Returns symbolic expression
```
"""
function Jacobi_c(α, κ::AbstractVector{<:Integer}, s::AbstractVector{<:Integer}, m, α1, α2)
    println("🔍 Jacobi_c CALLED | α: $α | κ: $κ | s: $s | m: $m | α1: $α1 | α2: $α2")
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition κ: must be non-increasing"))
    end
    if !parvalid(s)
        throw(ArgumentError("Invalid partition s: must be non-increasing"))
    end
    
    # Base case: if κ = s, return 1
    if κ == s
        println("  → Base case: κ == s, returning 1")
        return Sym(1)
    end
    
    result = Sym(0)
    
    # Loop from 1 to min(length(s)+1, length(κ))
    println("  → Looping i from 1 to $(min(length(s) + 1, length(κ)))")
    for i in 1:min(length(s) + 1, length(κ))
        # Create si by incrementing s at position i
        if i == length(s) + 1
            si = [s..., 1]
            println("    → i=$i: si = [s..., 1] = $si")
        else
            si = copy(s)
            si[i] = si[i] + 1
            println("    → i=$i: si = s with si[$i] incremented = $si")
        end
        
        # Check if si is already sorted in descending order (Maple: si = sort(si, '>'))
        si_sorted = sort(si, rev=true)
        println("    → si_sorted: $si_sorted | si == si_sorted: $(si == si_sorted) | partition_ge(κ, si): $(partition_ge(κ, si))")
        if si == si_sorted && partition_ge(κ, si)
            # Remove trailing zeros for validation
            si_clean = filter(x -> x > 0, si)
            println("    → si_clean: $si_clean")
            if parvalid(si_clean)
                println("    → Valid si_clean, computing recursive term")
                gbc_val = GBC(α, κ, si_clean)
                println("    → gbc_val: $gbc_val")
                gbc_cont_val = GBC_cont_explicit(α, si_clean, i)
                println("    → gbc_cont_val: $gbc_cont_val")
                jacobi_c_rec = Jacobi_c(α, κ, si_clean, m, α1, α2)
                println("    → jacobi_c_rec: $jacobi_c_rec")
                term = gbc_val * gbc_cont_val * jacobi_c_rec
                println("    → term: $term")
                result = result + term
                println("    → cumulative result: $result")
            end
        end
    end
    
    # Compute sums
    ks = sum(κ)
    ss = sum(s)
    println("  → ks (sum of κ): $ks | ss (sum of s): $ss")
    
    # Compute denominator: ((α1+α2+2*(m-1)/α+2)*(ks-ss)+rho(α,κ)-rho(α,s))
    rho_k = rho(α, κ)
    rho_s = rho(α, s)
    val = (α1 + α2 + 2*(m - 1)/α + 2) * (ks - ss) + rho_k - rho_s
    println("  → rho(α, κ): $rho_k | rho(α, s): $rho_s")
    println("  → denominator val: $val")
    
    if val == 0
        println("  → val == 0, returning 0")
        return Sym(0)
    end
    
    result = result / val
    result = simplify(result)
    println("Jacobi_c | α: $α | κ: $κ | s: $s | m: $m | α1: $α1 | α2: $α2 | result: $result")
    return result
end

