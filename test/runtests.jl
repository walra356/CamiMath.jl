# SPDX-License-Identifier: MIT

# Copyright (c) 2023 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

using CamiMath
using Test

@testset "CamiMath.jl" begin

    println("CamiMath.jl | 166 runtests | runtume 40s (estimated) | start")
    
    @test Type_IOP(1, 1) == Int64 
    @test Type_IOP(big(1), 1) == BigInt
    @test Type_IOP(2, 1) == BigInt
    @test Type_IOP(2, 1, "test(2)"; msg=false) == BigInt
    @test isforward(fwd)
    @test !isforward(bwd)
    @test isbackward(bwd)
    @test !isbackward(fwd)
    @test !isregular(rev)
    @test isregular(reg)
    @test isreversed(rev)
    @test !isreversed(reg)
    
    @test sup(-5 // 2) == "⁻⁵ᐟ²"
    @test 'D' * sup("superscript") == "Dˢᵘᵖᵉʳˢᶜʳⁱᵖᵗ"
    @test sub(-5 // 2) == "₋₅⸝₂"
    @test sub("-5/2") == "₋₅⸝₂"
    @test sub("e") == "ₑ"
    @test frac(-5 // 2) == "-⁵/₂"
    @test strRational(-5) == "-5"
    @test strRational(-5 // 2) == "-5/2"

    @test_throws DomainError bernoulliB(-1)
    @test eltype(bernoulliB(35)) == Rational{Int}
    @test eltype(bernoulliB(big(35))) == Rational{BigInt}
    @test bernoulliB(10; arr=true) == (1 // 1, -1 // 2, 1 // 6, 0 // 1, -1 // 30, 0 // 1, 1 // 42, 0 // 1, -1 // 30, 0 // 1, 5 // 66)
    @test sum([Rational{BigInt}(bernoulliB(n, msg=false)) for n = 1:90]) == 17080392099538483734383296956025848377395298523292305346008923087528282372539916092669722994190437011 // 3961456982724258461775089600226385
    @test bernoulliB(60; msg=false) == bernoulliB(60; msg=false, arr=true)[end]
    #...........................................................................
    @test divisor([1 // 1, 1 // 2, 1 // 3]) == 6
    @test numerators([1 // 1, 1 // 2, 1 // 3]) == [6, 3, 2]
    @test normalize_rationals([1 // 1, 1 // 2, 1 // 3]) == ([6, 3, 2], 6)
    #...........................................................................
    @test typeof(bigfactorial(5)) == Int
    typeof(bigfactorial(big(5))) == BigInt
    @test bigfactorial(21; msg=false) == 51090942171709440000
 #   @test bigconvert([[1 // 1, 1 // 2], [1 // 1, 1 // 2]]) == Vector{Vector{Rational{BigInt}}}
    @test convertToBig([[1 // 1, 1 // 2], [1 // 1, 1 // 2]]) == [[1 // 1, 1 // 2], [1 // 1, 1 // 2]]
    #...........................................................................
    @test_throws DomainError faulhaber_polynom(-1)
    @test eltype(faulhaber_polynom(10)) == Rational{Int}
    @test eltype(faulhaber_polynom(big(10))) == Rational{BigInt}
    @test sum([sum(faulhaber_polynom(p; msg=false)) for p = 1:90]) == 90 // 1
    @test sum([sum(faulhaber_polynom(big(p); msg=false)) for p = 1:90]) == 90 // 1
    #...........................................................................
    @test_throws DomainError faulhaber_polynomial(-1, 2)
    @test_throws DomainError faulhaber_polynomial(1, -1)
    @test faulhaber_polynomial(0, 2) == 0
    @test faulhaber_polynomial(1, 0) == 0
    @test faulhaber_polynomial(5, 30; msg=false) == 186552813930161650665
    @test faulhaber_polynomial(5, 37; msg=false) == 14556637744944425468330179
    @test faulhaber_polynomial(3, 6) == 276
    @test faulhaber_summation(3, 5) == 276

    @test_throws DomainError harmonicNumber(-1)
    @test sum([harmonicNumber(big(n)) for n = 1:46]) == 217436794888004994869 // 1345655451257488800
    @test sum(harmonicNumber(big(46); arr=true)) == 217436794888004994869 // 1345655451257488800
    @test sum([big(harmonicNumber(n; msg=false)) for n = 1:50]) == 10904958651492685640759 // 60765578514627386400
    @test sum(harmonicNumber(50; arr=true, msg=false)) == 10904958651492685640759 // 60765578514627386400

    @test_throws DomainError harmonicNumber(-1, 2)
    @test harmonicNumber(3, -5) == faulhaber_summation(3, 5)
    @test typeof(harmonicNumber(12, 3)) == Rational{Int}
    @test typeof(harmonicNumber(big(12), 3)) == Rational{BigInt}
    @test sum([big(harmonicNumber(n, 1; msg=false)) for n = 1:50]) == 10904958651492685640759 // 60765578514627386400
    @test sum([big(harmonicNumber(n, 2; msg=false)) for n = 1:30]) == 249434823919965027461489839 // 5424658191543895143840000
    @test sum([big(harmonicNumber(n, 3; msg=false)) for n = 1:20]) == 946050957591455286003049 // 40049466474634102886400
    @test sum([big(harmonicNumber(n, 4; msg=false)) for n = 1:15]) == 33970914744683670577361 // 2107930685520179520000
    @test sum([big(harmonicNumber(n, 5; msg=false)) for n = 1:12]) == 202913008605111524948813 // 16366888723117363200000
    @test sum([big(harmonicNumber(n, 6; msg=false)) for n = 1:10]) == 520072576088259310303 // 51219253009612800000
    @test sum([big(harmonicNumber(n, 7; msg=false)) for n = 1:10]) == 6501704870565012768992057 // 645362587921121280000000
    @test sum([big(harmonicNumber(n, 8; msg=false)) for n = 1:10]) == 16322500356140762610539037347 // 1626313721561225625600000000
    @test sum([big(harmonicNumber(n, 9; msg=false)) for n = 1:10]) == 41056936708131458470001164762001 // 4098310578334288576512000000000
    @test sum([big(harmonicNumber(n, 10; msg=false)) for n = 1:10]) == 20673934657221575836904008710237871 // 2065548531480481442562048000000000
    @test sum([big(harmonicNumber(n, 11; msg=false)) for n = 1:10]) == 260374709040010874103433717206249507497 // 26025911496654066176281804800000000000
    @test sum([big(harmonicNumber(n, 12; msg=false)) for n = 1:10]) == 655998094465816276746306450592922582177267 // 65585296971568246764230148096000000000000
    @test sum(harmonicNumber(6, 10; arr=true)) == 3630965833785900323 // 604661760000000000
    @test sum(harmonicNumber(4, 12; msg=false, arr=true)) == 35670966225905 // 8916100448256
    @test sum(harmonicNumber(10, 10; msg=false, arr=true)) == 20673934657221575836904008710237871 // 2065548531480481442562048000000000
    @test sum(harmonicNumber(10, 12; msg=false, arr=true)) == 655998094465816276746306450592922582177267 // 65585296971568246764230148096000000000000

    @test_throws DomainError fibonacci(-1)
    @test fibonacci(0) == 0
    @test typeof(fibonacci(92)) == Int
    @test typeof(fibonacci(99; msg=false)) == BigInt
    sum(fibonacci(92; arr=true)) == sum([fibonacci(i) for i = 1:92])
    @test sum(fibonacci(big(99); arr=true)) == sum([big(fibonacci(i; msg=false)) for i = 1:99])

    @test_throws DomainError canonical_partitions(5, -1)
    @test_throws DomainError canonical_partitions(5, 6)
    @test canonical_partitions(5) == [[1, 1, 1, 1, 1], [2, 2, 1], [3, 2], [4, 1], [5]]
    @test canonical_partitions(5, 0, reg) ==  [[5], [4, 1], [3, 2], [2, 2, 1], [1, 1, 1, 1, 1]]
    @test canonical_partitions(7, 3) == [3, 3, 1]
    @test integer_partitions(5) == [[1, 1, 1, 1, 1], [2, 2, 1], [2, 1, 1, 1], [3, 2], [3, 1, 1], [4, 1], [5]]
    @test integer_partitions(7, 4; transpose=true) == [[2, 2, 2, 1], [3, 2, 1, 1], [4, 1, 1, 1]]

    ftab = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];
    @test lagrange_polynom(ftab, 1:4, fwd) ≈ [0.0, 0.0, 1.0, 0.0]
    @test lagrange_polynom(ftab, 1:4, bwd) ≈ [9.0, 6.0, 1.0, 0.0]

    @test_throws DomainError laguerre_polynom(-1)
    @test typeof(laguerre_polynom(18)) == NTuple{19,Rational{Int}}
    @test typeof(laguerre_polynom(19, msg=false)) == Vector{Rational{BigInt}}
    @test sum(laguerre_polynom(20; msg=false)) == -21032925955607701 // 128047474114560000
    @test isapprox(sum(generalized_laguerre_polynom(20, 0.0; msg=false)), -big(21032925955607701) / big(128047474114560000); atol=1.0e-75)
    @test laguerre_polynom(8) == (1 // 1, -8 // 1, 14 // 1, -28 // 3, 35 // 12, -7 // 15, 7 // 180, -1 // 630, 1 // 40320)

    @test_throws DomainError generalized_laguerre_polynom(-1, 3)
    @test typeof(generalized_laguerre_polynom(18, 0; msg=false)) == NTuple{19,Rational{Int64}}
    @test typeof(generalized_laguerre_polynom(18, big(0); msg=false)) == Vector{Rational{BigInt}}
    @test typeof(generalized_laguerre_polynom(18, 0.0; msg=false)) == Vector{Float64}
    @test typeof(generalized_laguerre_polynom(18, big(0.0); msg=false)) == Vector{BigFloat}
    @test typeof(generalized_laguerre_polynom(19, 0; msg=false)) == Vector{Rational{BigInt}}
    @test typeof(generalized_laguerre_polynom(19, 0.0; msg=false)) == Vector{BigFloat}
    @test sum(generalized_laguerre_polynom(20, 0; msg=false)) == -21032925955607701 // 128047474114560000
    @test sum(generalized_laguerre_polynom(20, 0.0; msg=false)) ≈ -big(21032925955607701) / big(128047474114560000)
    @test generalized_laguerre_polynom(8, 3; msg=false) == [165 // 1, -330 // 1, 231 // 1, -77 // 1, 55 // 4, -11 // 8, 11 // 144, -11 // 5040, 1 // 40320]

    @test laguerreL(10, 5) == 254927 // 145152
    @test laguerreL(10, 15 // 3) == 254927 // 145152
    @test laguerreL(10, 5.0) == 1.7562761794536978

    @test generalized_laguerreL(10, 5, 5) == -425219 // 145152
    @test generalized_laguerreL(10, 5, 5.0) == -2.9294739307867648

    @test_throws DomainError hermite_polynom(-1)
    @test typeof(hermite_polynom(19)) == NTuple{20,Int}
    @test typeof(hermite_polynom(20, msg=false)) == Vector{BigInt}
    @test sum(hermite_polynom(20; msg=false)) == 1107214478336
    @test hermite_polynom(8) == (1680, 0, -13440, 0, 13440, 0, -3584, 0, 256)
    @test hermiteH(8, 5) == 52065680

    @test_throws DomainError pascal_triangle(-1)
    @test typeof(pascal_triangle(50)) == Vector{Int}
    @test pascal_triangle(0) == [1]
    @test pascal_triangle(5) == [1, 5, 10, 10, 5, 1]
    @test pascal_next([1, 4, 6, 4, 1]) == [1, 5, 10, 10, 5, 1]
    @test sum(sum.(pascal_triangle(5; arr=true, msg=false))) == 62
    @test sum(sum.(pascal_triangle(30; arr=true, msg=false))) == 2147483646

    @test permutations_unique_count([[1, 2], [1, 2, 3, 4, 5]], 2) == 120

    @test o = [pochhammer.([x for x = 0:-1:-p], p) for p = 0:5] == [[1], [0, -1], [0, 0, 2], [0, 0, 0, -6], [0, 0, 0, 0, 24], [0, 0, 0, 0, 0, -120]]
    @test o = [pochhammer.([x for x = 0:p], p) for p = 0:5] == [[1], [0, 1], [0, 2, 6], [0, 6, 24, 60], [0, 24, 120, 360, 840], [0, 120, 720, 2520, 6720, 15120]]
    @test pochhammer(big(-1) // big(50), 20) == -21605762356630090481082546653745369902321614221999 // 9536743164062500000000000000000000

    @test_throws DomainError polynomial((1, 1, 1, 1, 1), 1; deriv=-2)
    @test typeof(polynomial((1, 1, 1, 1, 1), 1)) == Int
    @test typeof(polynomial((1, 1, 1, 1, 1), big(1))) == BigInt
    @test typeof(polynomial((1, 1, 1, 1, 1), 1.0)) == Float64
    @test typeof(polynomial((1, 1, 1, 1, 1), big(1.0))) == BigFloat
    @test polynomial((1, 1, 1, 1, 1), 1; deriv=-1) == 137 // 60
    @test polynomial((1, 1, 1, 1, 1), 1; deriv=0) == 5
    @test polynomial((1, 1, 1, 1, 1), 1; deriv=1) == 10
    @test polynomial((1, 1, 1, 1, 1), 1; deriv=4) == 24
    @test polynomial((1, 1, 1, 1, 1), 1; deriv=5) == 0

    @test polynom_product((1.0, -1), (1, 2, -3)) == [1.0, 1.0, -5.0, 3.0]
    @test polynom_product((1, 2, -3), (1.0, -1)) == [1.0, 1.0, -5.0, 3.0]
    @test polynom_product((1 // 2, -1), (1, 2, -3)) == [1 // 2, 0 // 1, -7 // 2, 3 // 1]
    @test polynom_product((1, 1), (1, -1, 2)) == [1, 0, 1, 2]
    @test polynom_product_expansion((1, -1, 1), (1, 1, -1, 1, 1, 1), 4) == [1, 0, -1, 3, -1]
    @test polynom_product_expansion((1, -1, 1), (1, 1, -1, 1, 1, 1), 5) == [1, 0, -1, 3, -1, 1]
    @test polynom_product_expansion((1, -1, 1), (1, 1, -1, 1, 1, 1), 6) == [1, 0, -1, 3, -1, 1]
    @test polynom_product_expansion((1, 1, -1, 1, 1, 1), (1, -1, 1), 4) == [1, 0, -1, 3, -1]
    @test polynom_product_expansion((1, 1, -1, 1, 1, 1), (1, -1, 1), 5) == [1, 0, -1, 3, -1, 1]
    @test polynom_product_expansion((1, 1, -1, 1, 1, 1), (1, -1, 1), 6) == [1, 0, -1, 3, -1, 1]
    @test polynom_product_expansion((1, -1, 1), (1, 1, -1), 3) == [1, 0, -1, 2]
    @test polynom_product_expansion((1, -1, 1), (1, 1, -1), 4) == [1, 0, -1, 2, -1]
    @test polynom_product_expansion((1, -1, 1), (1, 1, -1), 5) == [1, 0, -1, 2, -1]
    @test polynom_power((1, 1, 1), 0) == [1]
    @test polynom_power((1, 1, 1), 1) == (1, 1, 1)
    @test polynom_power((1, 1, 1), 2) == [1, 2, 3, 2, 1]
    @test polynom_power((1, 1, 1), 3) == [1, 3, 6, 7, 6, 3, 1]

    @test texp(1 // 1, 0 // 1, 5) == 163 // 60
    @test texp(1, 0, 5) == 163 // 60
    @test texp(1.0, 0.0, 5) ≈ 2.7166666666666663
    @test typeof(texp(1 // 1, 0 // 1, 5)) == Rational{Int}
    @test typeof(texp(1, 0, 5)) == Rational{Int}
    @test typeof(texp(big(1), big(0), 5)) == Rational{BigInt}
    @test typeof(texp(1.0, 0.0, 5)) == Float64
    @test typeof(texp(big(1.0), big(0.0), 5)) == BigFloat

    @test log10_characteristic.([3, 30, 300]) == [0, 1, 2]
    @test log10_mantissa.([3, 30, 300]) ≈ [0.47712125471966244, 0.4771212547196624, 0.4771212547196626]

    @test istriangle(3, 4, 5) == true
    @test istriangle(1 // 2, 1, 1.5) == true
    @test triangle_coefficient(3, 4, 5) == 1 // 180180
    @test triangle_coefficient(1 // 2, 1, 1.5) == 1 // 12
    @test threeJsymbol(3, 0, 4, -1, 5, 1; msg=false) ≈ -0.10964174397241236
    @test CGC(3, 0, 4, -1, 5, -1; msg=false) ≈ -0.36364052611670256
    @test CGC(3, 0, 4, -1, 5, -1) ≈ -0.36364052611670256

end
