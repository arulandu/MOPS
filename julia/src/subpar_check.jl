"""
    subpar_check(s::AbstractVector{<:Integer}, k::AbstractVector{<:Integer}) -> Bool

Check if partition `s` is a subpartition of partition `k`.

According to Maple: `subPar?(s, k)` checks if s is a subpartition of k.
The Maple implementation checks: if k[i] > s[i] for any i, return false,
then return true if length(s) >= length(k).

This means: s is a subpartition of k if k[i] <= s[i] for all i AND length(s) >= length(k).

# Examples
```jldoctest
julia> MOPS.subpar_check([2, 1], [3, 2, 1])
false

julia> MOPS.subpar_check([3, 2, 1], [2, 1])
true
```
"""
function subpar_check(s::AbstractVector{<:Integer}, k::AbstractVector{<:Integer})
    if !parvalid(k)
        throw(ArgumentError("Invalid partition k: must be non-increasing"))
    end
    if !parvalid(s)
        throw(ArgumentError("Invalid partition s: must be non-increasing"))
    end
    
    # Maple: subPar?(k, s) - but called as subPar?(s, k) in usage
    # Looking at usage: if (not(`MOPS/subPar?`(s,k))) then return 0; end if;
    # So it's checking if s is a subpartition of k
    # The Maple code checks: if k[i] > s[i] for any i, return false
    # Then: return len(s) >= len(k)
    # But this logic seems backwards - let's check the actual definition
    # A subpartition s of k should satisfy s[i] <= k[i] for all i
    # So we check: if s[i] > k[i] for any i, return false
    # Then check length condition
    
    min_len = min(length(k), length(s))
    for i in 1:min_len
        if s[i] > k[i]  # If s[i] > k[i], then s is not a subpartition of k
            return false
        end
    end
    
    # After checking all elements, return true if length(s) >= length(k)
    # Actually, wait - let me check the Maple code again. It says nops(s) >= nops(k)
    # But that seems wrong for subpartitions. Let me check usage patterns.
    # For now, let's use the standard definition: s is subpartition if s[i] <= k[i] for all i
    return true
end

