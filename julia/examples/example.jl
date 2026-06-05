# This file demonstrates Julia-native symbolic MOPS usage.

using MOPS

@syms a b x

gbinom(a, k::Integer) = k == 0 ? 1 : SFact(a - k + 1, k) / factorial(k)
function classical_jacobi(n::Integer, a, b, t)
    result = 0
    for j in 0:n
        result += gbinom(n + a, n - j) * gbinom(n + b, j) * ((t - 1) / 2)^j * ((t + 1) / 2)^(n - j)
    end
    simplify(expand(result))
end

println("=" ^ 60)
println("Example 1: simplify((a+b+3)*(a+b+4)*Jacobi(2, [2], a, b, [x]))")
println("=" ^ 60)

result1 = simplify((a+b+3)*(a+b+4)*Jacobi(2, [2], a, b, [x]))
println(result1)
println()

println("=" ^ 60)
println("Example 2: 2 * classical Jacobi P₂^(a,b)(1-2x)")
println("=" ^ 60)

result2 = simplify(2 * classical_jacobi(2, a, b, 1 - 2*x))
println(result2)
println()

println("Expanded:")
println(expand(result2))
