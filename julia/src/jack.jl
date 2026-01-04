"""
    Jack(α, κ::AbstractVector{<:Integer}, args...) -> Any

Compute the Jack polynomial for partition `κ`.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition
- `ar1`: (optional) Number of variables, list of evaluation points, or normalization scale
- `ar2`: (optional) Normalization scale

# Examples
```jldoctest
julia> using SymPy
julia> @syms a x
julia> MOPS.Jack(a, [2], [x])
# Returns symbolic expression
```
"""
function Jack(α, κ::AbstractVector{<:Integer}, args...)
    if isempty(κ)
        return 1
    end
    
    # Parse arguments following Maple logic
    nd = 0  # n defined flag
    scd = 0  # scale defined flag
    n = nothing
    sc = nothing
    scale = 1
    
    if length(args) >= 2
        # nargs >= 2: n=ar1, sc=ar2 (or just sc if ar1 is symbol)
        ar1 = args[1]
        if length(args) >= 2
            ar2 = args[2]
        else
            ar2 = nothing
        end
        
        if ar1 isa Symbol || ar1 isa String
            # ar1 is normalization (nargs = 1 or 2 with ar1 as symbol)
            scd = 1
            sc = ar1 isa String ? Symbol(ar1) : ar1
            nd = 0
            if length(args) >= 2
                # If ar2 exists, it's n
                n = ar2
                nd = 1
            end
        else
            # ar1 is n
            nd = 1
            n = ar1
            if length(args) >= 2
                sc = ar2
                if sc isa String
                    sc = Symbol(sc)
                end
                scd = 1
            else
                sc = 0
            end
        end
    elseif length(args) >= 1
        # nargs = 1: ar1 could be normalization or n
        ar1 = args[1]
        if ar1 isa Symbol || ar1 isa String
            # ar1 is normalization
            scd = 1
            sc = ar1 isa String ? Symbol(ar1) : ar1
            nd = 0
        else
            # ar1 is n
            nd = 1
            n = ar1
            sc = 0
        end
    end
    
    kl = sum(κ)
    
    # Compute normalization scalar
    if scd == 1 && sc !== nothing && sc != 0
        sc_str = string(sc)
        if startswith(sc_str, "sym")
            if endswith(sc_str, "J")
                # Maple: return `MOPS/JtoC`(a,k,0) * J[op(k)]
                return JtoC(α, κ, 0) * jack_sym(κ)
            elseif endswith(sc_str, "P")
                # Maple: return `MOPS/PtoC`(a,k,0) * P[op(k)]
                return PtoC(α, κ, 0) * p_sym(κ)
            else
                # symC or other
                if n isa AbstractVector
                    # Maple: return `MOPS/eval`(a,k,n)
                    # eval = evalJack / JtoC for C normalization
                    return evalJack(α, κ, n) / JtoC(α, κ, 0)
                else
                    # Maple: return C[op(k)]
                    return c_sym(κ)
                end
            end
        elseif sc == :J || sc == "J"
            # Maple: scale := `MOPS/JtoC`(a,k,0)
            scale = JtoC(α, κ, 0)
        elseif sc == :P || sc == "P"
            # Maple: scale := `MOPS/PtoC`(a,k,0)
            scale = PtoC(α, κ, 0)
        end
    end
    
    # Compute jc
    m = if n === nothing
        0
    elseif n isa Integer
        n
    elseif n isa AbstractVector
        length(n)
    else
        0
    end
    jc = JtoC(α, κ, m)
    
    # If n is a list, evaluate at those points
    if n isa AbstractVector
        return simplify(scale / jc * evalJack(α, κ, n))
    end
    
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    if nd == 1 && n isa Integer && length(κ) > n
        return 0
    end
    
    result = 0
    
    # Faster algorithm for single-part partitions
    if length(κ) == 1
        for i in Par(kl)
            if nd == 1 && n isa Integer && length(i) > n
                continue
            end
            p1 = 1
            p2 = 1
            for j in 1:length(i)
                p1 = p1 * factorial(i[j])
                p2 = p2 * SFact(1/α + 1, i[j] - 1)
            end
            # In Maple: m[op(i)] represents monomial symmetric function
            coeff = simplify(scale / jc * factorial(kl) / p1 * α^(kl - length(i)) * p2)
            # Multiply coefficient by monomial symmetric function m[i]
            result = result + coeff * monomial_sym(i)
        end
        return result
    end
    
    # General case: use dominate or Par
    if length(κ) <= 3 && kl / length(κ) > 4
        f = 1
        s = Par(kl)
    else
        f = 0
        s = dominate(κ)
    end
    
    for i in s
        # Compute coefficients
        if partition_le(i, κ) && (f == 0 || (f == 1 && isdominate(κ, i)))
            if nd == 0
                coeff = simplify(scale * Jack_c(α, κ, i))
            else
                coeff = simplify(scale * Jack_c(α, κ, i, m))
            end
            # Multiply coefficient by monomial symmetric function m[i]
            result = result + coeff * monomial_sym(i)
        end
    end
    
    return result
end

