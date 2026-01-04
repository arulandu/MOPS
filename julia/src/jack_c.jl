"""
    Jack_c(α, k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer}, m::Integer) -> Any

Compute the Jack polynomial coefficient for partition `k` in the expansion
of the monomial symmetric function `m[k]` in terms of Jack polynomials.

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `k`: The partition for the Jack polynomial
- `l`: The partition for the monomial symmetric function
- `m`: Number of variables (optional, defaults to sum(k))

# Examples
```jldoctest
julia> using SymPy
julia> @syms a
julia> MOPS.Jack_c(a, [2], [1, 1])
# Returns coefficient
```
"""
function Jack_c(α, k::AbstractVector{<:Integer}, l::AbstractVector{<:Integer}, m::Integer = sum(k))
    if !parvalid(k)
        throw(ArgumentError("Invalid partition k: must be non-increasing"))
    end
    if !parvalid(l)
        throw(ArgumentError("Invalid partition l: must be non-increasing"))
    end
    
    if length(l) > m
        return 0
    end
    
    # Base case: if k == l
    if k == l
        result = 1
        for j in 1:length(k)
            for i in 1:k[j]
                result = result * Uhook(α, k, i, j)
            end
        end
        ks = sum(k)
        return (α^ks * factorial(ks)) / result
    end
    
    result = 0
    u = copy(l)
    
    for i in 1:(length(l) - 1)
        for t in 1:(k[1] - u[i])
            br = false
            u = copy(l)
            u[i] = u[i] + t
            ui = copy(u)
            
            for j in (i+1):length(l)
                if br
                    break
                end
                u = copy(ui)
                u[j] = u[j] - t
                
                if u[j] < 0
                    br = true
                    break
                end
                
                # Sort in descending order
                sort!(u, rev=true)
                
                # Remove trailing zeros
                while !isempty(u) && u[end] == 0
                    pop!(u)
                end
                
                # Check if u <= l or u > k
                if partition_le(u, l) || partition_gt(u, k)
                    br = true
                    break
                end
                
                # Recursive call
                result = result + (l[i] - l[j] + 2*t) * Jack_c(α, k, u, m)
            end
        end
    end
    
    # Compute rho values
    i_rho = rho(α, k)
    j_rho = rho(α, l)
    
    if i_rho == j_rho
        return 0
    end
    
    return simplify((2/α) / (i_rho - j_rho) * result)
end

