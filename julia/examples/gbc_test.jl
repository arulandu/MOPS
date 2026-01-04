using MOPS
using SymPy

"""
    test_GBC_cont(α; max_n=5, verbose=true)

Compare GBC_cont and GBC_cont_explicit on all partitions of size ≤ max_n
and all valid contiguous indices n.

Prints mismatches and returns true iff all tests pass.
"""
function test_GBC_cont(α; max_n=5, verbose=true)
    ok = true

    # helper: generate all partitions of k in nonincreasing form
    function partitions(k, maxpart=k)
        k == 0 && return [[]]
        res = Vector{Vector{Int}}()
        for p in min(k, maxpart):-1:1
            for tail in partitions(k-p, p)
                push!(res, [p; tail])
            end
        end
        return res
    end

    for k in 1:max_n
        for u in partitions(k)
            for n in 1:length(u)
                # contiguous move must be valid
                u[n] > 0 || continue

                # σ must be a partition after decrement
                σ = copy(u)
                σ[n] -= 1
                while !isempty(σ) && σ[end] == 0
                    pop!(σ)
                end
                parvalid(σ) || continue

                try
                    println("STARTING")
                    v1 = simplify(GBC_cont(α, u, n))
                    v2 = simplify(GBC_cont_explicit(α, u, n))
                    if simplify(v1 - v2) != 0
                        ok = false
                        println("❌ mismatch")
                        println("  u = ", u)
                        println("  n = ", n)
                        println("  GBC_cont         = ", v1)
                        println("  GBC_cont_explicit = ", v2)
                        println()
                        if verbose
                            return false
                        end
                    end
                catch e
                    ok = false
                    println("⚠️ error for u=$u, n=$n")
                    println(e)
                    if verbose
                        return false
                    end
                end
            end
        end
    end

    if ok
        println("✅ all tests passed up to size ", max_n)
    end
    return ok
end

@syms a
test_GBC_cont(a, max_n=4)
# GBC_cont_explicit(a, [2], 1)