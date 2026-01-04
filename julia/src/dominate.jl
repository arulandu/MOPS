"""
    dominate(μ::AbstractVector{<:Integer}) -> Vector{Vector{Int}}
    dominate(μ::AbstractVector{<:Integer}, n::Integer) -> Vector{Vector{Int}}

Generate all partitions dominated by partition `μ`.

If called with one argument, returns all partitions dominated by `μ`.
If called with two arguments, returns partitions dominated by `μ` with size at most `n`.

# Arguments
- `μ`: The partition
- `n`: (optional) Maximum size of partitions to generate

# Examples
```jldoctest
julia> MOPS.dominate([3, 1])
# Returns all partitions dominated by [3, 1]
```
"""
function dominate(μ::AbstractVector{<:Integer})
    if !parvalid(μ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    n = sum(μ)
    return dominate(μ, n)
end

function dominate(μ::AbstractVector{<:Integer}, n::Integer)
    if !parvalid(μ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    total = sum(μ)
    # In Maple: n = min(total, args[2]) if nargs > 1
    # But n is used as max number of parts: min(n, m+1) < j
    # However, for filtering by sum, we'll use n as max sum
    n_max_parts = min(total, n)  # For part count constraint
    n_max_sum = n  # For sum constraint (if n < total, filter by sum)
    
    res = [copy(μ)]
    sat = 1  # Index of next partition to process (1-based in Julia)
    
    while sat <= length(res)
        nu = res[sat]
        m = length(nu)
        nu0 = vcat(nu, [0])  # Append 0 for easier indexing
        
        for i in 1:m
            # Skip if consecutive elements are equal
            if nu0[i] == nu0[i+1]
                continue
            end
            
            # Determine j
            if nu0[i+1] + 1 < nu0[i]
                j = i + 1
            else
                # Find j such that nu0[j-1] == nu0[j]
                j = i + 2
                while j <= m && nu0[j-1] == nu0[j]
                    j += 1
                end
            end
            
            # Check if j is valid (constraint on number of parts)
            if min(n_max_parts, m + 1) < j
                continue
            end
            
            # Generate new partition
            if j <= m
                # subsop(i = nu[i]-1, j = nu[j]+1, nu)
                nu_new = copy(nu)
                nu_new[i] = nu[i] - 1
                nu_new[j] = nu[j] + 1
                # Ensure partition is still valid (non-increasing)
                sort!(nu_new, rev=true)
            else
                # [op(subsop(i = nu[i]-1, nu)), 1]
                nu_new = copy(nu)
                nu_new[i] = nu[i] - 1
                # Remove trailing zeros
                while !isempty(nu_new) && nu_new[end] == 0
                    pop!(nu_new)
                end
                push!(nu_new, 1)
                sort!(nu_new, rev=true)
            end
            
            # Check if new partition is valid and not already in res
            if parvalid(nu_new)
                # Check if not already in res
                if !any(ν -> ν == nu_new, res)
                    push!(res, nu_new)
                end
            end
        end
        
        sat += 1
    end
    
    # Filter by sum if n < total
    if n < total
        res = filter(ν -> sum(ν) <= n, res)
    end
    
    return res
end

