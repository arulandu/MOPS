@testset "Symbolics simplifier sanity" begin
    import SymbolicIntegration
    @syms a x

    # Direct Symbolics simplification should combine like rational terms.
    direct = Symbolics.simplify(2 / a + 4 / a)
    @test string(direct) == "6 / a"

    # Basic cancellation to exact zero through Symbolics directly.
    zero_direct = Symbolics.simplify((x + 1) - (x + 1))
    @test string(zero_direct) == "0"

    # Compare direct Symbolics simplify to the MOPS compatibility wrapper.
    wrapped = simplify(2 / a + 4 / a)
    @test string(direct) == string(wrapped)

    # Verify "weird zero-like literals": exact rational zero can be represented
    # as 0//1 in symbolic form (rather than plain "0").
    zero_rational = simplify(expand(Sym(0//1)))
    @test symbolic_zero(zero_rational)
    @test replace(string(Symbolics.value(Symbolics.Num(zero_rational))), " " => "") == "0//1"

    # Reproduce symbolic rational-zero forms seen in integration output and
    # ensure symbolic_zero handles them directly.
    @syms y
    function _sanity_poly_defint(expr, var, lo, hi)
        primitive = SymbolicIntegration.integrate(expand(expr), var)
        integral = Symbolics.substitute(primitive, Dict(var => hi)) - Symbolics.substitute(primitive, Dict(var => lo))
        return simplify(expand(integral))
    end
    poly = expand(Jacobi(1//1, [2], 0, 0, [x, y]))
    int_const = _sanity_poly_defint(_sanity_poly_defint(expand(poly * (x - y)^2), x, 0, 1), y, 0, 1)
    raw_int_const = MOPS._raw_symbolic(int_const)
    @test replace(string(raw_int_const), " " => "") in ("0//1", "(0//1)")
    @test symbolic_zero(int_const)

    # Verify near-zero floating constants are not symbolic zero.
    # This is the residual pattern seen when Float64 literals leak into algebra.
    near_zero = simplify(expand(0.1 + 0.2 - 0.3))
    @test !symbolic_zero(near_zero)
    @test isempty(Symbolics.get_variables(near_zero))
    @test abs(Float64(Symbolics.value(Symbolics.Num(near_zero)))) < 1e-12
end
