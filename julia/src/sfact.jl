"""
    SFact(a, b)

Compute the shifted factorial (Pochhammer symbol) (a)_b = a(a+1)(a+2)...(a+b-1).

This is equivalent to Γ(a+b)/Γ(a).

# Arguments
- `a`: Can be a number or a SymPy symbolic expression
- `b`: An integer (the number of factors)

# Examples
```jldoctest
julia> MOPS.SFact(3, 5)
5040

julia> using SymPy
julia> @vars n
julia> MOPS.SFact(n, 3)
n*(n + 1)*(n + 2)
```
"""
function SFact(a, b::Integer)
    if b < 0
        throw(DomainError(b, "b must be non-negative"))
    end
    
    if b == 0
        return 1
    end
    
    # Compute product: a * (a+1) * ... * (a+b-1)
    result = 1
    for i in 0:(b-1)
        result = result * (a + i)
    end
    return result
end
