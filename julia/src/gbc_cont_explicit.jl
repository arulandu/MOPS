"""
    GBC_cont_explicit(α, u::AbstractVector{<:Integer}, n::Integer) -> Any

Return the contiguous generalized binomial coefficient obtained by deleting one
box from row `n` of `u`.

Historically this Julia port had a second, hook-formula implementation here.
That implementation misinterpreted the column condition in the printed formula
and, for example, returned `3*(α+1)/2` for `(3 choose 2)` instead of the
ordinary one-row value `3`. To keep a single authoritative implementation,
this compatibility wrapper now delegates to `GBC_cont`, the hook-ratio
algorithm inherited from the original MOPS code and used by `GBC`.
"""
function GBC_cont_explicit(α, u::AbstractVector{<:Integer}, n::Integer)
    return GBC_cont(α, u, n)
end
