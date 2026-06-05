"""
    monomial_sym(κ::AbstractVector{<:Integer})

Create a formal Julia-native symbolic representation of the monomial symmetric
function ``m_κ``.  The symbol is intentionally atomic: MOPS manipulates these as
basis elements, not as expanded symmetric polynomials.
"""
function monomial_sym(κ::AbstractVector{<:Integer})
    isempty(κ) && return 1
    return _formal_symbol("m", κ)
end

function jack_sym(κ::AbstractVector{<:Integer})
    isempty(κ) && return 1
    return _formal_symbol("J", κ)
end

function p_sym(κ::AbstractVector{<:Integer})
    isempty(κ) && return 1
    return _formal_symbol("P", κ)
end

function c_sym(κ::AbstractVector{<:Integer})
    isempty(κ) && return 1
    return _formal_symbol("C", κ)
end
