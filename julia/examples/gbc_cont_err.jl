using MOPS
using SymPy

@syms a

u = [3]
n = 1

println("GBC_cont(a, [3], 1)")
println(GBC_cont(a, u, n))
println("GBC_cont_explicit(a, [3], 1)")
println(GBC_cont_explicit(a, u, n))