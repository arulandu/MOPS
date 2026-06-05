# Exact Gram-Schmidt verification of Appendix A.9 in arXiv:2510.04422.
# Independent of MOPS/DES. This computes P^8_[2,2,2],3,4(x I_3) from
# the defining Jacobi inner product and normalization P(0_3)=1.
#
# This file uses only Base/LinearAlgebra plus SymPy for final symbolic display.
# It should also work as a runtests.jl helper, though I could not run Julia in
# this sandbox.

using LinearAlgebra
using SymPy

const Exp3 = NTuple{3,Int}

function addexp(a::Exp3, b::Exp3)::Exp3
    return (a[1]+b[1], a[2]+b[2], a[3]+b[3])
end

function mul_poly(A::Dict{Exp3,BigInt}, B::Dict{Exp3,BigInt})
    C = Dict{Exp3,BigInt}()
    for (ea, ca) in A, (eb, cb) in B
        e = addexp(ea, eb)
        C[e] = get(C, e, big(0)) + ca*cb
    end
    filter!(kv -> kv[2] != 0, C)
    return C
end

function linear_power(i::Int, j::Int, power::Int)
    # 1-based i,j; expand (x_i - x_j)^power
    D = Dict{Exp3,BigInt}()
    for k in 0:power
        exp = [0,0,0]
        exp[i] = k
        exp[j] = power-k
        coeff = big(binomial(power,k)) * ((-1)^(power-k))
        e = (exp[1], exp[2], exp[3])
        D[e] = get(D, e, big(0)) + coeff
    end
    return D
end

function one_minus_power(i::Int, power::Int)
    D = Dict{Exp3,BigInt}()
    for k in 0:power
        exp = [0,0,0]
        exp[i] = k
        coeff = big(binomial(power,k)) * ((-1)^k)
        e = (exp[1], exp[2], exp[3])
        D[e] = get(D, e, big(0)) + coeff
    end
    return D
end

function partitions_of(d::Int, maxlen::Int)
    out = Vector{Vector{Int}}()
    function rec(rem::Int, maxpart::Int, pref::Vector{Int})
        if rem == 0
            push!(out, copy(pref))
        elseif length(pref) < maxlen
            for p in min(maxpart,rem):-1:1
                rec(rem-p, p, [pref; p])
            end
        end
    end
    rec(d,d,Int[])
    return out
end

partitions_upto(dmax::Int, maxlen::Int) = reduce(vcat, [partitions_of(d,maxlen) for d in 0:dmax])

function unique_perms(v::Vector{Int})
    out = Set{Exp3}()
    for a in v, b in v, c in v
        # Small n=3 version: include only permutations by comparing counts.
        p = [a,b,c]
        ok = true
        for z in unique(v)
            ok &= count(==(z), p) == count(==(z), v)
        end
        ok && push!(out, (a,b,c))
    end
    return collect(out)
end

function monomial_exponents(part::Vector{Int})
    padded = [part; zeros(Int, 3-length(part))]
    return unique_perms(padded)
end

function exact_a9_by_gram_schmidt()
    n = 3
    beta = 8
    gamma1 = 3
    gamma2 = 4
    target = [2,2,2]
    degree = sum(target)

    W = Dict{Exp3,BigInt}((0,0,0)=>big(1))
    for i in 1:n, j in i+1:n
        W = mul_poly(W, linear_power(i,j,beta))
    end
    for i in 1:n
        W = mul_poly(W, one_minus_power(i,gamma2))
    end
    println("expanded polynomial weight terms: ", length(W))

    moment_cache = Dict{Exp3,Rational{BigInt}}()
    function moment_monom(exp::Exp3)
        haskey(moment_cache, exp) && return moment_cache[exp]
        total = big(0)//big(1)
        for (mon, coeff) in W
            den = big(exp[1] + mon[1] + gamma1 + 1) *
                  big(exp[2] + mon[2] + gamma1 + 1) *
                  big(exp[3] + mon[3] + gamma1 + 1)
            total += coeff//den
        end
        moment_cache[exp] = total
        return total
    end

    mex_cache = Dict{Tuple{Int,Int,Int},Vector{Exp3}}()
    function exps_for(part::Vector{Int})
        keyv = [part; zeros(Int, 3-length(part))]
        key = (keyv[1],keyv[2],keyv[3])
        return get!(mex_cache, key) do
            monomial_exponents(part)
        end
    end

    inner_cache = Dict{Tuple{Tuple{Int,Int,Int},Tuple{Int,Int,Int}},Rational{BigInt}}()
    function keypart(p::Vector{Int})
        v = [p; zeros(Int, 3-length(p))]
        return (v[1],v[2],v[3])
    end
    function inner(p::Vector{Int}, q::Vector{Int})
        key = (keypart(p), keypart(q))
        haskey(inner_cache, key) && return inner_cache[key]
        total = big(0)//big(1)
        for a in exps_for(p), b in exps_for(q)
            total += moment_monom((a[1]+b[1], a[2]+b[2], a[3]+b[3]))
        end
        inner_cache[key] = total
        return total
    end

    lower = partitions_upto(degree-1, n)
    println("lower-degree basis size: ", length(lower))

    G = Matrix{Rational{BigInt}}(undef, length(lower), length(lower))
    rhs = Vector{Rational{BigInt}}(undef, length(lower))
    for i in eachindex(lower)
        for j in eachindex(lower)
            G[i,j] = inner(lower[j], lower[i])
        end
        rhs[i] = -inner(target, lower[i])
    end

    coeffs = G \ rhs

    @syms x
    function mult(part::Vector{Int})
        return length(exps_for(part))
    end
    diag = Sym(mult(target)) * x^degree
    for (a, mu) in zip(coeffs, lower)
        diag += SymPy.Sym(a) * mult(mu) * x^sum(mu)
    end
    P = simplify(expand(diag / subs(diag, x=>0)))

    A9_24738 = (7315*x^6 - 19950*x^5 + 24738*x^4 - 17472*x^3 + 7488*x^2 - 1872*x + 208) / Sym(208)
    A9_printed = (7315*x^6 - 19950*x^5 + 24762*x^4 - 17472*x^3 + 7488*x^2 - 1872*x + 208) / Sym(208)

    println("P = ", P)
    println("difference from 24738 version = ", simplify(expand(P - A9_24738)))
    println("difference from printed 24762 version = ", simplify(expand(P - A9_printed)))

    return P
end

exact_a9_by_gram_schmidt()
