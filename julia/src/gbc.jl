"""
    GBC(α, k::AbstractVector{<:Integer}, s::AbstractVector{<:Integer}) -> Any

Generalized binomial coefficients for Jack polynomials.

Computes the generalized binomial coefficient GBC(α, k, s) which appears
in the expansion of Jack polynomials.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `k`: The larger partition
- `s`: The smaller partition (subpartition of k)

# Examples
```jldoctest
julia> using SymPy
julia> @syms a
julia> MOPS.GBC(a, [2, 1], [1])
# Returns symbolic expression
```
"""
function GBC(α, k::AbstractVector{<:Integer}, s::AbstractVector{<:Integer})
    println("🔍 GBC CALLED | α: $α | k: $k | s: $s")
    if !parvalid(k)
        throw(ArgumentError("Invalid partition k: must be non-increasing"))
    end
    if !parvalid(s)
        throw(ArgumentError("Invalid partition s: must be non-increasing"))
    end
    
    # Base cases
    if k == s
        println("🔍 GBC: k == s case, returning 1")
        return Sym(1)
    end
    
    if isempty(s) || s == [0]
        println("🔍 GBC: isempty(s) || s == [0] case, returning 1")
        return Sym(1)
    end
    
    if s == [1]
        println("🔍 GBC: s == [1] case, returning sum(k)")
        return sum(k)
    end
    
    # Check if s is a subpartition of k
    if !subpar_check(s, k)
        println("🔍 GBC: !subpar_check(s, k) case, returning 0")
        return Sym(0)
    end
    
    result = Sym(0)
    
    # Loop from 1 to min(length(k), length(s)+1)
    for i in 1:min(length(k), length(s) + 1)
        # Create si by incrementing s at position i
        if i == length(s) + 1
            si = [s..., 1]
        else
            si = copy(s)
            si[i] = si[i] + 1
        end
        
        # Check if si == k (Maple checks this first)
        if si == k
            println("🔍 GBC: si == k case, calling GBC_cont_explicit | k: $k | i: $i")
            return simplify(GBC_cont_explicit(α, k, i))
        end
        
        # Check if si is already sorted in descending order (Maple: si = sort(si, '>'))
        si_sorted = sort(si, rev=true)
        if si == si_sorted && partition_ge(k, si)
            # Remove trailing zeros for validation
            si_clean = filter(x -> x > 0, si)
            if parvalid(si_clean)
                println("🔍 GBC: recursive case, calling GBC_cont_explicit | si_clean: $si_clean | i: $i")
                result = result + GBC(α, k, si_clean) * GBC_cont_explicit(α, si_clean, i)
            end
        end
    end
    
    # Compute sums
    ss = sum(s)
    ks = sum(k)
    
    if ks == ss
        return Sym(0)  # Avoid division by zero
    end
    
    result = result / (ks - ss)
    result = simplify(result)
    println("GBC | alpha: ", α, " | k: ", k, " | s: ", s, " | result: ", result)
    return result
end

