"""
    rho(α, κ::AbstractVector{<:Integer})

Compute the ρ function of a partition.

Given the Jack parameter `α` and the partition `κ`:
```math
ρ(α, κ) = \\sum_{i=1}^{\\operatorname{length}(κ)} κ_i \\left(κ_i - 1 - \\frac{2(i-1)}{α}\\right)
```

# Arguments
- `α`: The Jack parameter (can be numeric or symbolic)
- `κ`: The partition

# Examples
```jldoctest
julia> MOPS.rho(2, [3, 1]) |> Int
5
```
"""
function rho(α, κ::AbstractVector{<:Integer})
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    result = 0
    for i in 1:length(κ)
        result = result + κ[i] * (κ[i] - 1 - (2 / α) * (i - 1))
    end
    return result
end

