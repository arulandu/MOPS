"""
    subpar_check(s::AbstractVector{<:Integer}, k::AbstractVector{<:Integer}) -> Bool

Check if partition `s` is a subpartition (Young diagram inclusion) of partition `k`.

With partitions represented as row-lengths (non-increasing integer vectors),
`s ⊆ k` means

    sᵢ ≤ kᵢ  for all i ≥ 1,

where missing parts are treated as 0.

# Examples
```jldoctest
julia> MOPS.subpar_check([2, 1], [3, 2, 1])
true

julia> MOPS.subpar_check([3, 2, 1], [2, 1])
false
```
"""
function subpar_check(s::AbstractVector{<:Integer}, k::AbstractVector{<:Integer})
    if !parvalid(k)
        throw(ArgumentError("Invalid partition k: must be non-increasing"))
    end
    if !parvalid(s)
        throw(ArgumentError("Invalid partition s: must be non-increasing"))
    end
    
    # Check row-wise inclusion, padding missing parts of k with 0.
    for i in 1:length(s)
        ki = i <= length(k) ? k[i] : 0
        if s[i] > ki
            return false
        end
    end
    return true
end

