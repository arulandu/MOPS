"""
    egen(κ::AbstractVector{<:Integer}) -> Vector{Vector{Int}}

Generate valid eigenvalues (subpartitions) for evalJack.

This function generates a list of partitions that are used in the recursive
evaluation of Jack polynomials.

# Arguments
- `κ`: The partition (conjugate form)

# Examples
```jldoctest
julia> MOPS.egen([2, 1])
# Returns list of valid partitions
```
"""
function egen(κ::AbstractVector{<:Integer})
    if length(κ) == 1
        if κ[1] == 1
            return [Int[], [1]]
        else
            return [[κ[1] - 1], copy(κ)]
        end
    end
    
    # Recursive case: egen(k[1..nops(k)-1])
    lst1 = egen(κ[1:(length(κ)-1)])
    lst2 = Vector{Vector{Int}}()
    
    for l in lst1
        if length(l) == length(κ) - 1
            if l[end] < κ[end]
                l_new = vcat(l, [κ[end] - 1])
                if l_new[end] == 0
                    l_new = l_new[1:(end-1)]
                end
                push!(lst2, l_new)
            else
                l2 = vcat(l, [κ[end]])
                l_new = vcat(l, [κ[end] - 1])
                if l_new[end] == 0
                    l_new = l_new[1:(end-1)]
                end
                push!(lst2, l_new)
                push!(lst2, l2)
            end
        else
            push!(lst2, l)
        end
    end
    
    return lst2
end

