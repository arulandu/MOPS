function GBC_cont_explicit(α, u::AbstractVector{<:Integer}, n::Integer)
    # println("🔍 GBC_cont_explicit CALLED | α: $α | u: $u | n: $n")
    if !parvalid(u)
        throw(ArgumentError("Invalid partition u: must be non-increasing"))
    end
    if n < 1 || n > length(u)
        throw(ArgumentError("Index n must be between 1 and length(u)"))
    end
    if u[n] <= 0
        throw(ArgumentError("Cannot decrement u[n] since u[n] ≤ 0"))
    end

    # σ = u with n-th part decremented by 1 (trim trailing zeros)
    σ = collect(u)
    σ[n] -= 1
    while !isempty(σ) && σ[end] == 0
        pop!(σ)
    end
    if !parvalid(σ)
        throw(ArgumentError("Decrementing u at n=$(n) does not yield a valid partition σ"))
    end

    # IMPORTANT: paper's "i-th column" means i = n (the index incremented), not the added-box column
    distinguished_col = n
    # println("  → σ: $σ | distinguished_col: $distinguished_col")

    # Iterate Ferrers cells of κ in your coordinate convention: (a,b) = (col,row)
    cells(κ) = ((a,b) for b in 1:length(κ) for a in 1:κ[b])

    # j_σ = ∏ Uhook(σ) ∏ Lhook(σ)
    jσ = Sym(1)
    for (a,b) in cells(σ)
        jσ *= Uhook(α, σ, a, b)
        jσ *= Lhook(α, σ, a, b)
    end
    jσ = simplify(jσ)
    # println("  → jσ: $jσ")

    # g = (∏ A)(∏ B) with the MIXED hook choices from eq.(9)
    prodA = Sym(1)
    prodB = Sym(1)
    for (a,b) in cells(σ)
        if a == distinguished_col
            # in i-th column: A = LOWER(σ), B = UPPER(u)
            prodA *= Lhook(α, σ, a, b)
            prodB *= Uhook(α, u, a, b)
        else
            # not in i-th column: A = UPPER(σ), B = LOWER(u)
            prodA *= Uhook(α, σ, a, b)
            prodB *= Lhook(α, u, a, b)
        end
    end
    prodA = simplify(prodA)
    prodB = simplify(prodB)
    # println("  → prodA: $prodA")
    # println("  → prodB: $prodB")

    result = simplify((prodA * prodB) / jσ)
    # println("GBC_cont_explicit | α: $α | u: $u | n: $n | result: $result")
    return result
end
