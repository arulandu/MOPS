# Example from example.md
# This file demonstrates the Julia equivalents of the Maple examples

using MOPS
using SymPy

# Define symbolic variables
@syms a b x

println("=" ^ 60)
println("Example 1: simplify((a+b+3)*(a+b+4)*jacobi(2, [2], a, b, [x]))")
println("=" ^ 60)

# Compute the multivariate Jacobi polynomial
result1 = (a+b+3)*(a+b+4)*Jacobi(2, [2], a, b, [x])

println("Result:")
println(result1)
println()

println("Simplified:")
result1_simplified = simplify(result1)
println(result1_simplified)
println()

println("=" ^ 60)
println("Example 2: simplify(2 * JacobiP(2, a, b, (1-2x)))")
println("=" ^ 60)

# For the univariate Jacobi polynomial, we use SymPy's jacobi function
# Note: SymPy uses different parameter names - it's jacobi(n, alpha, beta, x)
# where alpha and beta are the parameters, and n is the degree
result2 = 2 * sympy.jacobi(2, a, b, 1 - 2*x)

println("Result:")
println(result2)
println()

println("Simplified:")
result2_simplified = simplify(result2)
println(result2_simplified)
println()

println("Grouped by x:")
result2_expanded = expand(result2_simplified)
result2_grouped = collect(result2_expanded, x)
println(result2_grouped)
println()

println("=" ^ 60)
println("Both examples completed successfully!")
println("=" ^ 60)
