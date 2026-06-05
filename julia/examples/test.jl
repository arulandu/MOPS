using MOPS

@syms a b x

res = Jacobi(2, [2], a, b, [x])
println(res)
