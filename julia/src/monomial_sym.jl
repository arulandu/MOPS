"""
    monomial_sym(κ::AbstractVector{<:Integer}) -> Sym

Create a symbolic representation of the monomial symmetric function m[κ].

In Maple, this is represented as `m[op(κ)]` where `op(κ)` expands the partition.

# Arguments
- `κ`: The partition

# Examples
```jldoctest
julia> using SymPy
julia> m = MOPS.monomial_sym([2])
# Returns m[2]

julia> m = MOPS.monomial_sym([1, 1])
# Returns m[1,1]
```
"""
# Global IndexedBase for monomial symmetric functions
const _m_base = Ref{Any}(nothing)

function _get_m_base()
    if _m_base[] === nothing
        _m_base[] = sympy.IndexedBase("m")
    end
    return _m_base[]
end

function monomial_sym(κ::AbstractVector{<:Integer})
    if isempty(κ)
        return 1
    end
    m_base = _get_m_base()
    # Create indexed symbol m[κ[1], κ[2], ...]
    # Use Python's __getitem__ method for indexing
    if length(κ) == 1
        return m_base.__getitem__(κ[1])
    else
        return m_base.__getitem__(Tuple(κ))
    end
end

# Also export a function to get J[κ], P[κ], C[κ] for symbolic representations
const _J_base = Ref{Any}(nothing)
const _P_base = Ref{Any}(nothing)
const _C_base = Ref{Any}(nothing)

function _get_J_base()
    if _J_base[] === nothing
        _J_base[] = sympy.IndexedBase("J")
    end
    return _J_base[]
end

function _get_P_base()
    if _P_base[] === nothing
        _P_base[] = sympy.IndexedBase("P")
    end
    return _P_base[]
end

function _get_C_base()
    if _C_base[] === nothing
        _C_base[] = sympy.IndexedBase("C")
    end
    return _C_base[]
end

function jack_sym(κ::AbstractVector{<:Integer})
    if isempty(κ)
        return 1
    end
    j_base = _get_J_base()
    if length(κ) == 1
        return j_base.__getitem__(κ[1])
    else
        return j_base.__getitem__(Tuple(κ))
    end
end

function p_sym(κ::AbstractVector{<:Integer})
    if isempty(κ)
        return 1
    end
    p_base = _get_P_base()
    if length(κ) == 1
        return p_base.__getitem__(κ[1])
    else
        return p_base.__getitem__(Tuple(κ))
    end
end

function c_sym(κ::AbstractVector{<:Integer})
    if isempty(κ)
        return 1
    end
    c_base = _get_C_base()
    if length(κ) == 1
        return c_base.__getitem__(κ[1])
    else
        return c_base.__getitem__(Tuple(κ))
    end
end

