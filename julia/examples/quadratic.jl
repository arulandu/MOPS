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

println("***Univariate***")
res = simplify((a+b+2)*Jacobi(2//1, [1], a, b, [x]))
expected = classical_jacobi(1, a, b, 1-2*x)
println("Result:")
println(res)
println("Expected:")
println(expected)
println("Result - Expected:")
println(simplify(expand(res - expected)))

println("***Quadratic***")
res = simplify((a+b+3)*(a+b+4)*Jacobi(2//1, [2], a, b, [x]))
expected = 2 * classical_jacobi(2, a, b, 1-2*x)
println("Result:")
println(res)
println("Expected:")
println(expected)
println("Result - Expected:")
println(simplify(expand(res - expected)))
