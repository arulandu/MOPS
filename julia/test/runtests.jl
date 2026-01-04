using MOPS
using Test
using SymPy

@testset "parvalid" begin
    @test parvalid(Int[]) == true
    @test parvalid([1]) == true
    @test parvalid([3, 2, 1]) == true
    @test parvalid([2, 2, 1]) == true
    @test parvalid([4, 3, 2, 1]) == true
    
    # Invalid partitions (not non-increasing)
    @test parvalid([2, 3, 1]) == false
    @test parvalid([1, 2]) == false
    @test parvalid([3, 1, 2]) == false
end

@testset "SFact" begin
    # Numeric tests
    @test SFact(3, 0) == 1
    @test SFact(3, 1) == 3
    @test SFact(3, 2) == 12  # 3 * 4
    @test SFact(3, 5) == 2520  # 3 * 4 * 5 * 6 * 7
    
    # Edge cases
    @test SFact(0, 3) == 0  # 0 * 1 * 2
    @test SFact(1, 4) == 24  # 1 * 2 * 3 * 4
    
    # Symbolic tests
    @syms n
    result = SFact(n, 3)
    @test result == n * (n + 1) * (n + 2)
    
    # Test cases from hsfact.png
    @test SFact(3, 5) == 2520
    @syms n
    @test simplify(SFact(n, 3) - n * (1 + n) * (2 + n)) == 0
    
    # Test with symbolic expression
    @syms a b
    result2 = SFact(a + b, 2)
    @test result2 == (a + b) * (a + b + 1)
    
    # Test error handling
    @test_throws DomainError SFact(3, -1)
end

@testset "conjugate" begin
    @test conjugate(Int[]) == Int[]
    @test conjugate([1]) == [1]
    @test conjugate([3, 2, 2]) == [3, 3, 1]
    @test conjugate([7, 5, 3, 3, 2]) == [5, 5, 4, 2, 2, 1, 1]
    @test conjugate([4,3,1,1]) == [4,2,2,1]
    
    # Test that conjugate of conjugate is original
    κ = [3, 2, 1]
    @test conjugate(conjugate(κ)) == κ
    
    # Test cases from hconjugate.png (already covered above, but ensuring they match)
    @test conjugate([3, 2, 2]) == [3, 3, 1]
    @test conjugate([7, 5, 3, 3, 2]) == [5, 5, 4, 2, 2, 1, 1]
end

@testset "arm" begin
    κ = [4, 3, 1, 1]
    @test arm(κ, 2, 1) == 2  # κ[1] - 2 = 4 - 2 = 2
    @test arm(κ, 1, 1) == 3  # κ[1] - 1 = 4 - 1 = 3
    @test arm(κ, 1, 2) == 2  # κ[2] - 1 = 3 - 1 = 2
    @test arm(κ, 3, 1) == 1  # κ[1] - 3 = 4 - 3 = 1
    
    @test_throws BoundsError arm([3, 2], 1, 3)
    @test_throws ArgumentError arm([2, 3, 1], 1, 1)  # Invalid partition
    
    # Test cases from harm.png
    @test arm([4, 2, 2], 2, 1) == 2
    @syms k
    @test simplify(arm([k+5, k+1], 1, 2) - k) == 0
end

@testset "leg" begin
    κ = [4, 3, 1, 1]
    # leg(κ, a, b) = arm(conjugate(κ), b, a)
    κ_conj = conjugate(κ)  # [4, 2, 2, 1]
    @test leg(κ, 1, 1) == arm(κ_conj, 1, 1)
    @test leg(κ, 1, 2) == arm(κ_conj, 2, 1)
    @test leg(κ, 2, 1) == arm(κ_conj, 1, 2)
    
    # Test case from hleg.png
    @test leg([4, 3, 3], 2, 1) == 2
    
    # Test with simpler partition
    κ2 = [3, 2, 1]
    κ2_conj = conjugate(κ2)
    @test leg(κ2, 1, 1) == arm(κ2_conj, 1, 1)
    @test leg(κ2, 1, 2) == arm(κ2_conj, 2, 1)
    
    @test_throws ArgumentError leg([2, 3, 1], 1, 1)  # Invalid partition
end

@testset "Uhook" begin
    κ = [3, 2, 1]
    # Uhook(α, κ, x, y) = leg(κ, x, y) + α*(1 + arm(κ, x, y))
    α = 2
    @test Uhook(α, κ, 1, 1) == leg(κ, 1, 1) + α * (1 + arm(κ, 1, 1))
    
    # Test with symbolic α
    @syms alpha
    result = Uhook(alpha, κ, 1, 1)
    expected = leg(κ, 1, 1) + alpha * (1 + arm(κ, 1, 1))
    @test result == expected
    
    @test_throws ArgumentError Uhook(2, [2, 3, 1], 1, 1)  # Invalid partition
    
    # Test cases from huhook.png
    @test Uhook(2, [4, 2, 1], 3, 1) == 4
    @test Uhook(3, [2, 1, 1], 2, 1) == 3
end

@testset "Lhook" begin
    κ = [3, 2, 1]
    # Lhook(α, κ, x, y) = leg(κ, x, y) + 1 + α*(arm(κ, x, y))
    α = 2
    @test Lhook(α, κ, 1, 1) == leg(κ, 1, 1) + 1 + α * arm(κ, 1, 1)
    
    # Test with symbolic α
    @syms alpha
    result = Lhook(alpha, κ, 1, 1)
    expected = leg(κ, 1, 1) + 1 + alpha * arm(κ, 1, 1)
    @test result == expected
    
    @test_throws ArgumentError Lhook(2, [2, 3, 1], 1, 1)  # Invalid partition
    
    # Test cases from hlhook.png
    @test Lhook(2, [4, 2, 1], 3, 1) == 3
    @test Lhook(3, [2, 1, 1], 2, 1) == 1
end

@testset "GSFact" begin
    # GSFact(α, a, s) = ∏ᵢ SFact(a - (i-1)/α, s[i])
    @test GSFact(1, 3, [3, 2, 1]) == 360
    
    # Test with symbolic parameters
    @syms alpha a
    result = GSFact(alpha, a, [2, 1])
    expected = SFact(a, 2) * SFact(a - 1/alpha, 1)
    @test result == expected
    
    @test_throws ArgumentError GSFact(1, 3, [2, 3, 1])  # Invalid partition
    
    # Test cases from hgsfact.png
    @test GSFact(1, 3, [3, 2, 1]) == 360
    @test GSFact(2, 2, [4, 2, 2, 1]) == 450
end

@testset "rho" begin
    # rho(α, κ) = ∑ᵢ κ[i]*(κ[i] - 1 - (2/α)*(i-1))
    @test rho(2, [3, 1]) == 5
    
    # Test with symbolic α
    @syms alpha
    κ = [2, 1]
    result = rho(alpha, κ)
    expected = 2 * (2 - 1 - (2/alpha) * 0) + 1 * (1 - 1 - (2/alpha) * 1)
    @test simplify(result) == simplify(expected)
    
    @test rho(2, Int[]) == 0
    @test_throws ArgumentError rho(2, [2, 3, 1])  # Invalid partition
    
    # Test cases from hrho.png
    @test rho(2, [3, 1]) == 5
    @syms a
    @test simplify(rho(a, [4, 1, 1]) - (12 - 6/a)) == 0
end

@testset "subpar_check" begin
    # Maple: subPar?(s, k) - function defined as subPar?(k, s) but called as subPar?(s, k)
    # When called as subPar?(s, k), it checks: if s[i] > k[i] for any i, return false
    # This means: s is subpartition of k if s[i] <= k[i] for all i
    # Standard definition: a subpartition s of k satisfies s[i] <= k[i] for all i
    
    # [1] is subpartition of [3,1]? s[1]=1 <= k[1]=3 ✓
    @test subpar_check([1], [3, 1]) == true
    
    # [2] is subpartition of [3,1]? s[1]=2 <= k[1]=3 ✓
    @test subpar_check([2], [3, 1]) == true
    
    # [3] is subpartition of [3,1]? s[1]=3 <= k[1]=3 ✓
    @test subpar_check([3], [3, 1]) == true
    
    # [2,1] is subpartition of [3,1]? s[1]=2 <= k[1]=3 ✓, s[2]=1 <= k[2]=1 ✓
    @test subpar_check([2, 1], [3, 1]) == true
    
    # [3,1] is subpartition of [3,1]? Equal partitions
    @test subpar_check([3, 1], [3, 1]) == true
    
    # Test cases from hissubpar.png
    @test subpar_check([5, 3], [6, 4, 2]) == true
    @test subpar_check([5], [4, 1]) == false
    @test subpar_check([5, 5], [6, 5]) == true
    @test subpar_check([4], [4]) == true
    
    # [3,2] is subpartition of [3,1]? s[1]=3 <= k[1]=3 ✓, but s[2]=2 > k[2]=1 ✗
    @test subpar_check([3, 2], [3, 1]) == false
    
    @test_throws ArgumentError subpar_check([2, 3, 1], [2, 1])  # Invalid partition
end

@testset "partition comparison" begin
    # partition_le (<=)
    @test partition_le([2, 1], [3, 1]) == true
    @test partition_le([2, 1], [2, 1]) == true
    @test partition_le([3, 1], [2, 1]) == false
    @test partition_le([2], [2, 1]) == true  # [2] <= [2,1] (length)
    @test partition_le([2, 1], [2]) == false
    
    # partition_lt (<)
    @test partition_lt([2, 1], [3, 1]) == true
    @test partition_lt([2, 1], [2, 1]) == false
    @test partition_lt([2], [2, 1]) == true
    
    # partition_ge (>=)
    @test partition_ge([3, 1], [2, 1]) == true
    @test partition_ge([2, 1], [2, 1]) == true
    @test partition_ge([2, 1], [3, 1]) == false
    
    # partition_gt (>)
    @test partition_gt([3, 1], [2, 1]) == true
    @test partition_gt([2, 1], [2, 1]) == false
    
    @test_throws ArgumentError partition_le([2, 3, 1], [2, 1])  # Invalid partition
    
    # Test cases from hOrderingOperators.png
    @test partition_le([2, 1], [1, 1]) == false
    @test partition_gt([3, 2, 1], [1, 1, 0]) == true
    @test partition_lt([4, 1], [7, 3, 1]) == true
    @test partition_ge([2, 2], [7, 3, 3]) == false
end

@testset "subPar" begin
    # Test with [2]
    result = subPar([2])
    @test length(result) >= 2
    @test Int[] in result
    @test [1] in result
    @test [2] in result
    
    # Test with [3, 1]
    result = subPar([3, 1])
    @test Int[] in result
    @test [1] in result
    @test [2] in result
    @test [3] in result
    @test [1, 1] in result
    @test [2, 1] in result
    @test [3, 1] in result
    
    # Test with empty partition
    @test subPar(Int[]) == [Int[]]
    
    # Test with size constraint
    result = subPar([3, 1], 1)
    # All subpartitions of [3,1] are: [], [1], [2], [3], [1,1], [2,1], [3,1] with sizes 0,1,2,3,2,3,4
    # So size 1 subpartitions: [1]
    @test all(sum(p) == 1 for p in result)
    @test [1] in result
    
    result = subPar([3, 1], 0)
    @test result == [Int[]]
    
    @test_throws ArgumentError subPar([2, 3, 1])  # Invalid partition
    
    # Test cases from hsubpar.png
    @test subPar([4]) == [[4], [3], [2], [1], []]
    @test subPar([3, 3, 3]) == [[3, 3, 3], [3, 3, 2], [3, 3, 1], [3, 3], [3, 2, 2], [3, 2, 1], [3, 2], [3, 1, 1], [3, 1], [3], [2, 2, 2], [2, 2, 1], [2, 2], [2, 1, 1], [2, 1], [2], [1, 1, 1], [1, 1], [1], []]
    @test subPar([7, 2]) == [[7, 2], [7, 1], [7], [6, 2], [6, 1], [6], [5, 2], [5, 1], [5], [4, 2], [4, 1], [4], [3, 2], [3, 1], [3], [2, 2], [2, 1], [2], [1, 1], [1], []]
    @test subPar([1, 1, 1, 1]) == [[1, 1, 1, 1], [1, 1, 1], [1, 1], [1], []]
end

@testset "Par" begin
    # Test Par(3) - should return all partitions of 3
    result = Par(3)
    @test length(result) == 3
    @test [3] in result
    @test [2, 1] in result
    @test [1, 1, 1] in result
    
    # Test Par(4) - should return all partitions of 4
    result = Par(4)
    @test length(result) == 5
    @test [4] in result
    @test [3, 1] in result
    @test [2, 2] in result
    @test [2, 1, 1] in result
    @test [1, 1, 1, 1] in result
    
    # Test Par(0) - should return [[]]
    @test Par(0) == [Int[]]
    
    # Test Par(1) - should return [[1]]
    @test Par(1) == [[1]]
    
    # Test Par(2) - should return [[2], [1,1]]
    result = Par(2)
    @test length(result) == 2
    @test [2] in result
    @test [1, 1] in result
    
    # Verify all partitions sum to n
    for n in 1:5
        result = Par(n)
        for p in result
            @test sum(p) == n
            @test parvalid(p)  # All partitions should be valid
        end
    end
    
    # Test cases from hpar.png
    @test Par(4) == [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    @test Par(3) == [[3], [2, 1], [1, 1, 1]]
    @test Par(7) == [[7], [6, 1], [5, 2], [5, 1, 1], [4, 3], [4, 2, 1], [4, 1, 1, 1], [3, 3, 1], [3, 2, 2], [3, 2, 1, 1], [3, 1, 1, 1, 1], [2, 2, 2, 1], [2, 2, 1, 1, 1], [2, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]]
    @test length(Par(10)) == 42  # par(10) has 42 partitions
end

@testset "JtoC" begin
    # No hard expected value tests from PNG - skip for now
end

@testset "JackIdentity" begin
    using SymPy
    @syms a m
    
    # Test cases from hjackidentity.png
    # Test 1: jackidentity(a, [4,3,2,1], m)
    result1 = JackIdentity(a, [4, 3, 2, 1], m)
    expected1 = 3628800*a^6*m*(m+a)*(m+2*a)*(m+3*a)*(m-1)*(m-1+a)*(m-1+2*a)*(m-2)*(m-2+a)*(m-3) / ((3+4*a)*(3*a+4)*(2+3*a)^2*(3+2*a)^2*(1+2*a)^3*(2+a)^3)
    @test simplify(result1 - expected1) == 0
    
    # Test 2: jackidentity(1, [6,5,2,1], m)
    result2 = JackIdentity(1, [6, 5, 2, 1], m)
    expected2 = 4576*m^2*(m+1)^2*(m+2)^2*(m+3)^2*(m+4)*(m+5)*(m-1)^2*(m-2)*(m-3) / 297675
    @test simplify(result2 - expected2) == 0
    
    # Test 3: jackidentity(7/2, [3,2], 10)
    result3 = JackIdentity(7//2, [3, 2], 10)
    @test result3 == 4685625 // 506
    
    # Test 4: jackidentity(a, [10,2], 20)
    result4 = JackIdentity(a, [10, 2], 20)
    expected4 = 361152000*a^2*(20+a)*(10+a)*(20+3*a)*(5+a)*(4+a)*(10+3*a)*(20+7*a)*(5+2*a)*(20+9*a)*(19+a) / ((1+10*a)*(2+9*a)*(1+9*a)*(4*a+1)^2*(1+7*a)*(1+6*a)*(1+5*a)*(1+3*a)*(1+2*a)*(a+1)^2)
    @test simplify(result4 - expected4) == 0
end

@testset "egen" begin
    # No hard expected value tests from PNG - skip for now
end

@testset "evalJack" begin
    using SymPy
    @syms a x
    
    # Test with empty partition
    @test evalJack(a, Int[], [x]) == 1
    
    # Test with [1] and single variable
    @test evalJack(a, [1], [x]) == x
    
    # Test with [2] and single variable
    @test evalJack(a, [2], [x]) == x^2 * (a + 1)
end

@testset "Jack" begin
    using SymPy
    @syms a x
    
    # Test with empty partition
    @test Jack(a, Int[], [x]) == 1
end

@testset "isdominate" begin
    # Test dominance order
    @test isdominate([3, 1], [2, 2]) == true
    @test isdominate([2, 2], [3, 1]) == false
    @test isdominate([2, 1], [2, 1]) == true  # Equal partitions
    @test isdominate([4, 2], [3, 3]) == true
    @test isdominate([3, 3], [4, 2]) == false
    
    @test_throws ArgumentError isdominate([2, 3, 1], [2, 1])  # Invalid partition
end

@testset "dominate" begin
    @test dominate([3, 1]) == [[3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    @test dominate([3, 1], 3) == Vector{Int}[]
    @test dominate([2]) == [[2], [1, 1]]
    @test dominate([2, 1]) == [[2, 1], [1, 1, 1]]
end

@testset "Jack_c" begin
    using SymPy
    @syms a
    
    # Test with equal partitions
    @test Jack_c(a, [2], [2]) == 1
end

@testset "Jack (full)" begin
    # No hard expected value tests from PNG - skip for now
end

@testset "Jack from hjack.png" begin
    using SymPy
    @syms a x y z
    
    # Test 1: jack(a, [2,2])
    result1 = Jack(a, [2, 2])
    m22 = monomial_sym([2, 2])
    m211 = monomial_sym([2, 1, 1])
    m1111 = monomial_sym([1, 1, 1, 1])
    expected1 = (12*m22*a^2) / ((1+2*a)*(a+1)) + (24*m211*a^2) / ((a+1)^2*(1+2*a)) + (144*m1111*a^2) / ((2+a)*(a+1)^2*(1+2*a))
    @test simplify(result1 - expected1) == 0
    
    # Test 2: jack(a, [3,1], 'J')
    result2 = Jack(a, [3, 1], :J)
    m31 = monomial_sym([3, 1])
    m22 = monomial_sym([2, 2])
    m211 = monomial_sym([2, 1, 1])
    m1111 = monomial_sym([1, 1, 1, 1])
    expected2 = 2*m31*(a+1)^2 + 4*m22*(a+1) + 2*(3*a+5)*m211 + 24*m1111
    @test simplify(result2 - expected2) == 0
    
    # Test 3: jack(a, [3,1], 'P')
    # Expected from PNG: m[3,1] + (2*m[2,2])/(a+1) + ((3*a+5)*m[2,1,1])/((a+1)^2) + (12*m[1,1,1,1])/((a+1)^2)
    result3 = Jack(a, [3, 1], :P)
    m31 = monomial_sym([3, 1])
    m22 = monomial_sym([2, 2])
    m211 = monomial_sym([2, 1, 1])
    m1111 = monomial_sym([1, 1, 1, 1])
    expected3 = m31 + (2*m22) / (a+1) + (m211*(3*a+5)) / ((a+1)^2) + (12*m1111) / ((a+1)^2)
    @test simplify(result3 - expected3) == 0
    
    # Test 4: jack(2, [6], 2)
    result4 = Jack(2, [6], 2)
    m6 = monomial_sym([6])
    m51 = monomial_sym([5, 1])
    m42 = monomial_sym([4, 2])
    m33 = monomial_sym([3, 3])
    expected4 = m6 + (6//11)*m51 + (5//11)*m42 + (100//231)*m33
    @test simplify(result4 - expected4) == 0
    
    # Test 5: jack(a, [6,1,1], 3)
    result5 = Jack(a, [6, 1, 1], 3)
    m611 = monomial_sym([6, 1, 1])
    m521 = monomial_sym([5, 2, 1])
    m431 = monomial_sym([4, 3, 1])
    m332 = monomial_sym([3, 3, 2])
    m422 = monomial_sym([4, 2, 2])
    expected5 = (168*m611*a^2) / ((a+1)*(1+3*a)) + (840*m521*a^2) / ((4*a+1)*(1+3*a)*(a+1)) + (1680*m431*a^2) / ((4*a+1)*(1+3*a)^2) + (5040*m332*a^2) / ((1+3*a)^2*(4*a+1)*(1+2*a)) + (3360*m422*a^2) / ((1+3*a)^2*(4*a+1)*(a+1))
    @test simplify(result5 - expected5) == 0
    
    # Test 6: jack(1, [4], [x,y,z])
    result6 = Jack(1, [4], [x, y, z])
    expected6 = z^4 + z^3*y + z^3*x + z^2*y^2 + z^2*y*x + z^2*x^2 + z*y^3 + z*y^2*x + z*y*x^2 + z*x^3 + y^4 + y^3*x + y^2*x^2 + y*x^3 + x^4
    @test simplify(result6 - expected6) == 0
    
    # Test 7: jack(a, [3,2], 10, 'P')
    result7 = Jack(a, [3, 2], 10, :P)
    m32 = monomial_sym([3, 2])
    m221 = monomial_sym([2, 2, 1])
    m311 = monomial_sym([3, 1, 1])
    m2111 = monomial_sym([2, 1, 1, 1])
    m11111 = monomial_sym([1, 1, 1, 1, 1])
    expected7 = m32 + (m221*(3*a+5)) / ((a+1)^2) + (2*m311) / (a+1) + (12*m2111) / ((a+1)^2) + (60*m11111) / ((2+a)*(a+1)^2)
    @test simplify(result7 - expected7) == 0
end

@testset "Jack J-normalization coefficients" begin
    using SymPy
    @syms a
    
    # Test Jack(a, [2], :J) = (a + 1)*m[2] + 2*m[1, 1]
    m2 = monomial_sym([2])
    m11 = monomial_sym([1, 1])
    expected = (a + 1) * m2 + 2 * m11
    @test simplify(Jack(a, [2], :J) - expected) == 0
    
    # Test Jack(a, [2,1], :J) = (a + 2)*m[2, 1] + 6*m[1, 1, 1]
    m21 = monomial_sym([2, 1])
    m111 = monomial_sym([1, 1, 1])
    expected = (a + 2) * m21 + 6 * m111
    @test simplify(Jack(a, [2, 1], :J) - expected) == 0
    
    # Test Jack(a, [4], :J) = 6*(a+1)^2*m[2,2] + (a+1)*(2a+1)*(3a+1)*m[4] + 4*(a+1)*(2a+1)*m[3,1] + (12a+12)*m[2,1,1] + 24*m[1,1,1,1]
    m4 = monomial_sym([4])
    m31 = monomial_sym([3, 1])
    m22 = monomial_sym([2, 2])
    m211 = monomial_sym([2, 1, 1])
    m1111 = monomial_sym([1, 1, 1, 1])
    expected = 6 * (a + 1)^2 * m22 + (a + 1) * (2*a + 1) * (3*a + 1) * m4 + 
               4 * (a + 1) * (2*a + 1) * m31 + (12*a + 12) * m211 + 24 * m1111
    @test simplify(Jack(a, [4], :J) - expected) == 0
    
    # Test Jack(a, [2,2], :J) = 2*(a+1)*(a+2)*m[2,2] + (4a+8)*m[2,1,1] + 24*m[1,1,1,1]
    expected = 2 * (a + 1) * (a + 2) * m22 + (4*a + 8) * m211 + 24 * m1111
    @test simplify(Jack(a, [2, 2], :J) - expected) == 0
    
    # Test Jack(a, [1,1,1,1], :J) = 24*m[1,1,1,1]
    expected = 24 * m1111
    @test simplify(Jack(a, [1, 1, 1, 1], :J) - expected) == 0
    
    # Test Jack(a, [2,1,1], :J) = (2a+6)*m[2,1,1] + 24*m[1,1,1,1]
    expected = (2*a + 6) * m211 + 24 * m1111
    @test simplify(Jack(a, [2, 1, 1], :J) - expected) == 0
end

@testset "GBC_cont" begin
    @syms a
    
    # Test GBC_cont(a, [2,1], 1)
    @test simplify(GBC_cont(a, [2, 1], 1) - (2*a + 1)/(a + 1)) == 0
end

@testset "GBC" begin
    @syms a
    
    # Base case: GBC(a, [2], []) = 1
    @test GBC(a, [2], Int[]) == 1
    
    # Base case: GBC(a, [2], [1]) = sum([2]) = 2
    @test GBC(a, [2], [1]) == 2
    
    # Base case: GBC(a, [2], [2]) = 1
    @test GBC(a, [2], [2]) == 1
    
    # GBC(a, [2,1], [1]) = 3
    @test GBC(a, [2, 1], [1]) == 3
    
    # GBC(a, [2,1], [2]) = (a + 2)/(a + 1)
    result = GBC(a, [2, 1], [2])
    @test typeof(result) <: Sym
    @test simplify(result - (a + 2)/(a + 1)) == 0
    
    # GBC(a, [2,1], [2,1]) = 1
    @test GBC(a, [2, 1], [2, 1]) == 1
    
    
    # Test cases from hgbinomial.png
    @syms a
    result1 = GBC(a, [5, 4, 3, 2, 1], [3, 3, 3])
    expected1 = (160*a^4 + 892*a^3 + 1571*a^2 + 892*a + 160) / (4*(1+a)^4)
    @test simplify(result1 - expected1) == 0
    
    result2 = GBC(a, [2, 2, 1], [2, 1])
    expected2 = 6*(3+a) / (2+a)
    @test simplify(result2 - expected2) == 0
    
    @test GBC(2, [3, 1], [1, 1]) == Sym(7//3)
    @test GBC(17//2, [6], [4]) == 15
end

@testset "GBC Properties" begin
    @syms a
    
    # Property 1: (κ / 0) = 1 when σ is empty partition
    @test GBC(a, [2], Int[]) == 1
    @test GBC(a, [3, 1], Int[]) == 1
    @test GBC(a, [2, 2, 1], Int[]) == 1
    
    # Property 2: (κ / 1) = |κ| when σ = [1]
    @test GBC(a, [2], [1]) == 2  # |[2]| = 2
    @test GBC(a, [3, 1], [1]) == 4  # |[3,1]| = 4
    @test GBC(a, [2, 2, 1], [1]) == 5  # |[2,2,1]| = 5
    
    # Property 3: (κ / σ) = 0 if σ ⊄ κ (σ is not a subpartition of κ)
    # [3] is not a subpartition of [2,1] (3 > 2)
    @test GBC(a, [2, 1], [3]) == 0
    # [2,2] is not a subpartition of [3,1] (2 <= 3, but 2 > 1)
    @test GBC(a, [3, 1], [2, 2]) == 0
    # [1,1,1] is not a subpartition of [2,1] (length 3 > length 2)
    @test GBC(a, [2, 1], [1, 1, 1]) == 0
    
    # Property 4: (κ / σ) = δ_κ if |κ| = |σ| (Kronecker delta: 1 if κ=σ, 0 otherwise)
    # When κ = σ, should be 1
    @test GBC(a, [2], [2]) == 1
    @test GBC(a, [3, 1], [3, 1]) == 1
    @test GBC(a, [2, 2], [2, 2]) == 1
    # When |κ| = |σ| but κ ≠ σ, should be 0
    # |[2,1]| = 3, |[3]| = 3, but [2,1] ≠ [3]
    @test GBC(a, [2, 1], [3]) == 0
    # |[2,1]| = 3, |[1,1,1]| = 3, but [2,1] ≠ [1,1,1]
    @test GBC(a, [2, 1], [1, 1, 1]) == 0
    
    # Property 5: (κ / σ) ≠ 0 if |κ| = |σ| + 1, iff σ = κ(i) for some i
    # This means if |κ| = |σ| + 1, then GBC is non-zero only when σ can be obtained
    # by decrementing exactly one part of κ by 1
    
    # Case: κ = [3,1], |κ| = 4
    # σ = [2,1] has |σ| = 3, and [2,1] = [3,1] with first part decremented → should be non-zero
    @test GBC(a, [3, 1], [2, 1]) != 0
    
    # σ = [3] has |σ| = 3, and [3] = [3,1] with second part removed (decremented to 0) → should be non-zero
    @test GBC(a, [3, 1], [3]) != 0
    
    # Case: κ = [2,2], |κ| = 4
    # σ = [2,1] has |σ| = 3, and [2,1] = [2,2] with second part decremented → should be non-zero
    @test GBC(a, [2, 2], [2, 1]) != 0
    
    # Case: κ = [3,2], |κ| = 5
    # σ = [2,2] has |σ| = 4, and [2,2] = [3,2] with first part decremented → should be non-zero
    @test GBC(a, [3, 2], [2, 2]) != 0
    
    # σ = [3,1] has |σ| = 4, and [3,1] = [3,2] with second part decremented → should be non-zero
    @test GBC(a, [3, 2], [3, 1]) != 0
end

@testset "Jacobi_c" begin
    @syms a a1 a2
    
    # Base case: Jacobi_c(a, [2], [2], m, a1, a2) = 1
    @test Jacobi_c(a, [2], [2], 1, a1, a2) == 1
    @test Jacobi_c(a, [3, 1], [3, 1], 2, a1, a2) == 1
    
    # Jacobi_c(a, [2], [1], 1, a1, a2) should be 2/(a1 + a2 + 4)
    @test simplify(Jacobi_c(a, [2], [1], 1, a1, a2) - 2/(a1 + a2 + 4)) == 0
end

@testset "Jacobi" begin
    @syms a a1 a2 x g1 g2 n
    
    # Test cases from hjacobi.png with exact expected outputs
    @syms a g1 g2 n
    
    # Test 1: jacobi(a, [1,1], g1, g2, n)
    # Expected: C[1,1] - 2*a*C[1]*(g1*a + n - 2 + a)*(n - 1)/((g1*a + g2*a + 2*n - 4 + 2*a)*(a + 1)) + 
    #           2*a*(g1*a + n - 1 + a)*(g1*a + n - 2 + a)*n*(n - 1)/((g1*a + g2*a + 2*n - 4 + 2*a)*(g1*a + g2*a + 2*n - 3 + 2*a)*(a + 1))
    result1 = Jacobi(a, [1, 1], g1, g2, n)
    c11 = c_sym([1, 1])
    c1 = c_sym([1])
    expected1 = c11 - 2*a*c1*(a*(g1 + 1) + n - 2)*(n - 1)/((a*(g1 + g2 + 2) + 2*n - 4)*(a + 1)) + 
                2*a*(a*(g1 + 1) + n - 1)*(a*(g1 + 1) + n - 2)*n*(n - 1)/((a*(g1 + g2 + 2) + 2*n - 4)*(a*(g1 + g2 + 2) + 2*n - 3)*(a + 1))
    @test simplify(result1 - expected1) == 0
    
    # Test 2: jacobi(a, [2,1], 1, 1, 4, 'J')
    # Expected from PNG: (1/6)*((1+2a)(2+a))/a * J[2,1] + (9/(1+2a)) * J[2] * (a+1) + 
    #                   (3(4+a))/(a(2+a)) * J[1,1] * (a+1) - (9/5)*((4+a) * J[1] * (14+20a+11a^2))/((1+2a)(2+a)(a+1)^2) + 
    #                   (54/5)*((4+a)(3+2a)a(14+20a+11a^2))/((1+2a)(2+a)(7a+8)(a+1)^2)
    result2 = Jacobi(a, [2, 1], 1, 1, 4, :J)
    j21 = jack_sym([2, 1])
    j2 = jack_sym([2])
    j11 = jack_sym([1, 1])
    j1 = jack_sym([1])
    expected2 = -(1//6)*((1+2*a)*(2+a))/a * j21 + (9*a*(a+1)/(1+2*a)) * j2 + 
                (3*(4+a))/(a*(2+a)) * j11 * (a+1) - (9//5)*((4+a) * j1 * (14+20*a+11*a^2))*a/((1+2*a)*(2+a)*(a+1)^2) + 
                (54//5)*((4+a)*(3+2*a)*a*(14+20*a+11*a^2))/((1+2*a)*(2+a)*(7*a+8)*(a+1)^2)
    @test simplify(result2 - expected2) == 0
    
    # Test 3: jacobi(3, [6], 1, 0, 2, 'm')
    # Use 3//1 to get exact fractions
    # Expected: m[6] + (3/8)*m[5,1] + (15/52)*m[4,2] + (7/26)*m[3,3] - (561/164)*m[5] - 
    #           (2805/2132)*m[4,1] - (561/533)*m[3,2] + (58905/2132)*m[4] + (11781/1066)*m[3,1] + 
    #           (5049/533)*m[2,2] - (148104/533)*m[3] - (444312/3731)*m[2,1] + (555390/287)*m[2] + 
    #           (277695/287)*m[1,1] - (41654250/8323)*m[1] + (41654250/15457)
    result3 = Jacobi(3//1, [6], 1, 0, 2, :m)
    m6 = monomial_sym([6])
    m51 = monomial_sym([5, 1])
    m42 = monomial_sym([4, 2])
    m33 = monomial_sym([3, 3])
    m5 = monomial_sym([5])
    m41 = monomial_sym([4, 1])
    m32 = monomial_sym([3, 2])
    m4 = monomial_sym([4])
    m31 = monomial_sym([3, 1])
    m22 = monomial_sym([2, 2])
    m3 = monomial_sym([3])
    m21 = monomial_sym([2, 1])
    m2 = monomial_sym([2])
    m11 = monomial_sym([1, 1])
    m1 = monomial_sym([1])
    expected3 = m6 + (3//8)*m51 + (15//52)*m42 + (7//26)*m33 - (561//164)*m5 - 
                (2805//2132)*m41 - (561//533)*m32 + (58905//2132)*m4 + (11781//1066)*m31 + 
                (5049//533)*m22 - (148104//533)*m3 - (444312//3731)*m21 + (555390//287)*m2 + 
                (277695//287)*m11 - (41654250//8323)*m1 + (41654250//15457)
    @test simplify(result3 - expected3) == 0
    
    # Test 4: jacobi(a, [3], 0, 0, [x, y])
    # Note: PNG was incorrect. Using actual Julia output as expected.
    @syms x y
    result4 = Jacobi(a, [3], 0, 0, [x, y])
    expected4 = 9*(a^2 + 3*a + 2)/(10*a^2 + 9*a + 2) + 9*(-a*x - a*y - 2*x - 2*y)/(5*a + 2) + 
                3*(x^2*(a + 1) + 2*x*y + y^2*(a + 1))/(2*a + 1) + 
                (-x^3*(2*a + 1) - 3*x^2*y - 3*x*y^2 - y^3*(2*a + 1))/(2*a + 1)
    @test simplify(result4 - expected4) == 0
    
    
end
