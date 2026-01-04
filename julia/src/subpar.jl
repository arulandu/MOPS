"""
    subPar(κ::AbstractVector{<:Integer}) -> Vector{Vector{Int}}
    subPar(κ::AbstractVector{<:Integer}, n::Integer) -> Vector{Vector{Int}}

Generate all subpartitions of a partition.

If called with one argument, returns all subpartitions of `κ`.
If called with two arguments, returns all subpartitions of `κ` with size `n`.

# Examples
```jldoctest
julia> MOPS.subPar([2])
2-element Vector{Vector{Int64}}:
 []
 [1]

julia> MOPS.subPar([3, 1])
4-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3, 1]
```
"""
function subPar(κ::AbstractVector{<:Integer})
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    if isempty(κ)
        return [Int[]]
    end
    
    l = length(κ)
    
    # Find first position i where κ[i] != κ[i+1] (or i = l-1 if all equal)
    i = 1
    while i < l && κ[i] == κ[i+1]
        i += 1
    end
    
    # Determine j: if κ[i] > 1, j = κ[i] - 1, else nothing
    j = κ[i] > 1 ? κ[i] - 1 : nothing
    
    result = Vector{Vector{Int}}()
    
    # Case 1: Keep first i elements as κ[1], recurse on rest [op(i+1..l, mu)]
    # Maple: [seq([`$`(mu[1],i), op(nu)],nu = `MOPS/subPar`([op(i+1 .. l,mu)]))]
    if i < l
        rest = κ[(i+1):l]
        rest_subpars = subPar(rest)
        for nu in rest_subpars
            prefix = fill(κ[1], i)
            push!(result, vcat(prefix, nu))
        end
        # If rest is empty or subPar(rest) contains [], we also need [] in the final result
        # (Maple's structure includes [] from the recursive calls)
        if isempty(rest) || Int[] in rest_subpars
            # [] will be added at the end during deduplication
        end
    else
        # If i == l, all elements are the same
        # [op(i+1..l, mu)] = [op(l+1..l, mu)] = [] (empty range)
        # subPar([]) = [[]], so we get [[κ[1] repeated i times, op([])]] = [[κ[1] repeated i times]]
        rest = Int[]  # Empty list
        for nu in subPar(rest)  # subPar([]) = [[]]
            prefix = fill(κ[1], i)
            push!(result, vcat(prefix, nu))
        end
        # Only generate subpartitions with fewer elements when j = NULL (can't decrease)
        # Otherwise, smaller partitions will come from Case 2's recursion
        # For example: [2,2] with j=1 will get [2] from subPar([2,1]) in Case 2
        # But [1,1] with j=NULL needs [1] from Case 1's special handling
        if j === nothing && i > 1
            # Generate subpartitions by taking fewer elements
            for num_elements in (i-1):-1:1
                prefix = fill(κ[1], num_elements)
                push!(result, prefix)
            end
        end
    end
    
    # Case 2: Decrease element at position i: subsop(i = j, mu)
    # Maple: op(`MOPS/subPar`(subsop(i = j,mu)))
    if j !== nothing
        mu_new = copy(κ)
        mu_new[i] = j
        # Remove trailing zeros if any (to keep partition valid)
        while !isempty(mu_new) && mu_new[end] == 0
            pop!(mu_new)
        end
        append!(result, subPar(mu_new))
    end
    
    # Remove duplicates while preserving order (matching Maple's behavior)
    # The empty partition [] should be included if any recursive call would produce it
    # This happens when we can recursively reduce the partition to []
    unique_result = Vector{Vector{Int}}()
    seen = Set{Vector{Int}}()
    
    for r in result
        if !(r in seen)
            push!(unique_result, r)
            push!(seen, r)
        end
    end
    
    # Add empty partition at the end (Maple includes it from the deepest recursive calls)
    # The empty partition represents taking no elements, which is always a valid subpartition
    if !(Int[] in unique_result)
        push!(unique_result, Int[])
    end
    
    return unique_result
end

function subPar(κ::AbstractVector{<:Integer}, n::Integer)
    if !parvalid(κ)
        throw(ArgumentError("Invalid partition: must be non-increasing"))
    end
    
    m = sum(κ)
    
    if m < n || n < 0
        return Vector{Vector{Int}}()
    end
    
    if n == 0
        return [Int[]]  # Only the empty partition has size 0
    end
    
    if n == m
        return [copy(κ)]
    end
    
    l = length(κ)
    
    # Find first position where κ[i] != κ[i+1]
    i = 1
    while i < l && κ[i] == κ[i+1]
        i += 1
    end
    
    # Determine j: if κ[i] > 1, j = κ[i] - 1, else nothing
    j = κ[i] > 1 ? κ[i] - 1 : nothing
    
    result = Vector{Vector{Int}}()
    
    # Case 1: Keep first i elements as κ[1], recurse on rest
    # Maple: [seq([`$`(mu[1],i), op(nu)],nu = `MOPS/subPar`([op(i+1 .. l,mu)],n-i*mu[1]))]
    if i < l
        rest = κ[(i+1):l]
        for nu in subPar(rest, n - i * κ[1])
            prefix = fill(κ[1], i)
            push!(result, vcat(prefix, nu))
        end
    else
        # If i == l, all elements are the same
        # [op(i+1..l, mu)] is empty, so subPar([], n - i*mu[1])
        # If n - i*mu[1] == 0, we get [[]], otherwise []
        if n - i * κ[1] == 0
            prefix = fill(κ[1], i)
            push!(result, prefix)
        end
        # Also generate subpartitions with fewer elements when all are equal
        # For [1,1] with n=1, we should get [1] by taking 1 element
        if i > 1
            for num_elements in (i-1):-1:1
                if num_elements * κ[1] == n
                    prefix = fill(κ[1], num_elements)
                    push!(result, prefix)
                end
            end
        end
    end
    
    # Case 2: Decrease element at position i
    # Maple: op(`MOPS/subPar`(subsop(i = j,mu),n))
    if j !== nothing
        mu_new = copy(κ)
        mu_new[i] = j
        # Remove trailing zeros if any
        while !isempty(mu_new) && mu_new[end] == 0
            pop!(mu_new)
        end
        append!(result, subPar(mu_new, n))
    end
    
    return result
end

