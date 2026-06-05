# Julia-native symbolic compatibility layer for the MOPS port.
#
# The original project used SymPy.jl for three pieces of functionality:
#   * constructing scalar symbolic variables in examples/tests,
#   * simplifying/expanding rational symbolic expressions, and
#   * creating formal indexed symbols such as m[2,1].
# This file supplies the small API surface MOPS needs using Symbolics.jl and
# SymbolicUtils.jl only; it intentionally does not depend on Python.

import Symbolics
import SymbolicUtils

const Sym = Symbolics.Num

macro syms(args...)
    return esc(Expr(:macrocall, GlobalRef(Symbolics, Symbol("@variables")), __source__, args...))
end

function simplify(x; kwargs...)
    (x isa Number && !(x isa Symbolics.Num)) && return x
    try
        return Symbolics.simplify(x; kwargs...)
    catch
        return x
    end
end

function expand(x)
    (x isa Number && !(x isa Symbolics.Num)) && return x
    return Symbolics.expand(x)
end

function _raw_symbolic(x)
    x isa Symbolics.Num && return Symbolics.value(x)
    return x
end

function symbolic_zero(x)
    y = simplify(expand(x))
    raw = _raw_symbolic(y)
    if raw isa Number
        return iszero(raw)
    end
    # Constant symbolic expressions can appear as internal symbolic nodes
    # (e.g. "(0//1)") instead of plain Julia numbers.
    if isempty(Symbolics.get_variables(y))
        v = try
            Symbolics.value(Symbolics.Num(raw))
        catch
            return false
        end
        return v isa Number && iszero(v)
    end
    return false
end

function _formal_symbol(prefix::AbstractString, κ::AbstractVector{<:Integer})
    name = isempty(κ) ? Symbol(prefix) : Symbol(prefix * "[" * join(κ, ",") * "]")
    return Symbolics.variable(name)
end

function _symconvert(x)
    return x isa Symbolics.Num ? x : Symbolics.Num(x)
end
