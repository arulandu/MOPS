using MOPS
using SymPy

@syms a b c x y z g1 g2

# # univariate
# println("***Univariate***")
# res = simplify((a+b+2)*Jacobi(2//1, [1], a, b, [x]))
# expected = sympy.jacobi(1, a, b, 1-2*x)
# expected = collect(expand(simplify(expected)), x)
# println("Result:")
# println(res)
# println("Expected:")
# println(expected)
# println("Result - Expected:")
# println(simplify(res - expected))

# multivariate
println("***Multivariate***")
res = simplify((a+b+3)*(a+b+4)*Jacobi(2//1, [2], a, b, [x]))
expected = 2*sympy.jacobi(2, a, b, 1-2*x)
expected = collect(expand(simplify(expected)), x)
println("Result:")
println(res)
println("Expected:")
println(expected)
println("Result - Expected:")
println(simplify(res - expected))