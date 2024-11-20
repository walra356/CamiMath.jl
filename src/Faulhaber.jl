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

# ==============================================================================
#                          Faulhaber.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#               faulhaber_polynom(p::Integer [; msg=true])
# ------------------------------------------------------------------------------

function _faulhaber_BigInt(p::Integer)

    p = big(p)

    P = CamiMath.pascal_triangle(p)
    B = CamiMath.bernoulliB(p; arr=true)

    F = (B .* P) // p

    F[2] = -F[2]                      # equivalent to B[2] = -B[2]

    F = Base.reverse(F)               # reverse to standard order

    F[1] = big(0) // big(1)           # set polynomial constant: c_0 = 0//1


    return Tuple(F)

end

# ..............................................................................
@doc raw"""
    faulhaber_polynom(p::Integer [; msg=true])

Vector representation of the coefficients of the [`faulhaber_polynomial`](@ref) 
of degree `p`, 
```math
   c=[c_0,⋯\ c_p],
```
where ``c_0=0,\ \ c_j=\frac{1}{p}{\binom{p}{p-j}}B_{p-j}``,
with ``j∈\{ 1,⋯\ p\}``. The ``B_{p-j}`` are [`bernoulliB`](@ref) in the 
*even index convention* (but with 
``B_1=+\frac{1}{2}`` rather than ``-\frac{1}{2}``).

- `msg` : integer-overflow protection (IOP) - warning on activation 
(for `p > 36`)
#### Example:
```
faulhaber_polynom(6)
7-element Vector{Rational{Int64}}:
  0//1
  0//1
 -1//12
  0//1
  5//12
  1//2
  1//6
```
"""
function faulhaber_polynom(p::Integer; msg=true)

    N = (
        (0, 1), (0, 1, 1), (0, 1, 3, 2), (0, 0, 1, 2, 1),
        (0, -1, 0, 10, 15, 6), (0, 0, -1, 0, 5, 6, 2),
        (0, 1, 0, -7, 0, 21, 21, 6), (0, 0, 2, 0, -7, 0, 14, 12, 3),
        (0, -3, 0, 20, 0, -42, 0, 60, 45, 10),
        (0, 0, -3, 0, 10, 0, -14, 0, 15, 10, 2),
        (0, 5, 0, -33, 0, 66, 0, -66, 0, 55, 33, 6),
        (0, 0, 10, 0, -33, 0, 44, 0, -33, 0, 22, 12, 2),
        (0, -691, 0, 4550, 0, -9009, 0, 8580, 0, -5005, 0, 2730, 1365, 210),
        (0, 0, -691, 0, 2275, 0, -3003, 0, 2145, 0, -1001, 0, 455, 210, 30),
        (0, 105, 0, -691, 0, 1365, 0, -1287, 0, 715, 0, -273, 0, 105, 45, 6),
        (0, 0, 420, 0, -1382, 0, 1820, 0, -1287, 0, 572, 0, -182,
            0, 60, 24, 3),
        (0, -3617, 0, 23800, 0, -46988, 0, 44200, 0, -24310, 0, 8840,
            0, -2380, 0, 680, 255, 30),
        (0, 0, -10851, 0, 35700, 0, -46988, 0, 33150, 0, -14586, 0, 4420,
            0, -1020, 0, 255, 90, 10),
        (0, 219335, 0, -1443183, 0, 2848860, 0, -2678316, 0, 1469650,
            0, -529074, 0, 135660, 0, -27132, 0, 5985, 1995, 210),
        (0, 0, 438670, 0, -1443183, 0, 1899240, 0, -1339158, 0, 587860,
            0, -176358, 0, 38760, 0, -6783, 0, 1330, 420, 42),
        (0, -3666831, 0, 24126850, 0, -47625039, 0, 44767800, 0, -24551230,
            0, 8817900, 0, -2238390, 0, 426360, 0, -65835, 0, 11550, 3465, 330),
        (0, 0, -3666831, 0, 12063425, 0, -15875013, 0, 11191950, 0, -4910246,
            0, 1469650, 0, -319770, 0, 53295, 0, -7315, 0, 1155, 330, 30),
        (0, 4272565, 0, -28112371, 0, 55491755, 0, -52160757, 0, 28601650,
            0, -10266878, 0, 2600150, 0, -490314, 0, 72105, 0, -8855,
            0, 1265, 345, 30),
        (0, 0, 51270780, 0, -168674226, 0, 221967020, 0, -156482271,
            0, 68643960, 0, -20533756, 0, 4457400, 0, -735471, 0, 96140,
            0, -10626, 0, 1380, 360, 30), (0, -1181820455, 0, 7776068300,
            0, -15349354566, 0, 14427856300, 0, -7911048145, 0, 2839363800,
            0, -718681460, 0, 135207800, 0, -19684665, 0, 2302300, 0, -230230,
            0, 27300, 6825, 546),
        (0, 0, -1181820455, 0, 3888034150, 0, -5116451522, 0, 3606964075,
            0, -1582209629, 0, 473227300, 0, -102668780, 0, 16900975,
            0, -2187185, 0, 230230, 0, -20930, 0, 2275, 546, 42),
        (0, 538845489, 0, -3545461365, 0, 6998461470, 0, -6578294814,
            0, 3606964075, 0, -1294535151, 0, 327618900, 0, -61601268,
            0, 8947575, 0, -1036035, 0, 98670, 0, -8190, 0, 819, 189, 14),
        (0, 0, 1077690978, 0, -3545461365, 0, 4665640980, 0, -3289147407,
            0, 1442785630, 0, -431511717, 0, 93605400, 0, -15400317,
            0, 1988350, 0, -207207, 0, 17940, 0, -1365, 0, 126, 28, 2),
        (0, -23749461029, 0, 156265191810, 0, -308455138755, 0, 289936260900,
            0, -158975458005, 0, 57055613550, 0, -14439045915, 0, 2714556600,
            0, -394066935, 0, 45522750, 0, -4292145, 0, 339300, 0, -23751,
            0, 2030, 435, 30),
        (0, 0, -23749461029, 0, 78132595905, 0, -102818379585, 0, 72484065225,
            0, -31795091601, 0, 9509268925, 0, -2062720845, 0, 339319575,
            0, -43785215, 0, 4552275, 0, -390195, 0, 28275, 0, -1827,
            0, 145, 30, 2),
        (0, 8615841276005, 0, -56689963476223, 0, 111901503855141,
            0, -105183202315455, 0, 57673154564025, 0, -20698604632251,
            0, 5238144213225, 0, -984742931403, 0, 142933380975,
            0, -16502417085, 0, 1552325775, 0, -121486365, 0, 8099091,
            0, -484561, 0, 35805, 7161, 462),
        (0, 0, 68926730208040, 0, -226759853904892, 0, 298404010280376,
            0, -210366404630910, 0, 92277047302440, 0, -27598139509668,
            0, 5986450529400, 0, -984742931403, 0, 127051894200,
            0, -13201933668, 0, 1128964200, 0, -80990910, 0, 4984056,
            0, -276892, 0, 19096, 3696, 231),
        (0, -1780853160521127, 0, 11717544135366800, 0, -23129505098298984,
            0, 21740863606141680, 0, -11920762929084900, 0, 4278299465840400,
            0, -1082696242302360, 0, 203539317999600, 0, -29542287942090,
            0, 3410340318000, 0, -320618389080, 0, 25033554000, 0, -1652214564,
            0, 94143280, 0, -4869480, 0, 314160, 58905, 3570),
        (0, 0, -1780853160521127, 0, 5858772067683400, 0, -7709835032766328,
            0, 5435215901535420, 0, -2384152585816980, 0, 713049910973400,
            0, -154670891757480, 0, 25442414749950, 0, -3282476438010,
            0, 341034031800, 0, -29147126280, 0, 2086129500, 0, -127093428,
            0, 6724520, 0, -324632, 0, 19635, 3570, 210),
        (0, 90219075042845, 0, -593617720173709, 0, 1171754413536680,
            0, -1101405004680904, 0, 603912877948380, 0, -216741144165180,
            0, 54849993151800, 0, -10311392783832, 0, 1496612632350,
            0, -172761917790, 0, 16239715800, 0, -1267266360, 0, 83445180,
            0, -4707164, 0, 231880, 0, -10472, 0, 595, 105, 6),
        (0, 0, 541314450257070, 0, -1780853160521127, 0, 2343508827073360,
            0, -1652107507021356, 0, 724695453538056, 0, -216741144165180,
            0, 47014279844400, 0, -7733544587874, 0, 997741754900,
            0, -103657150674, 0, 8858026800, 0, -633633180, 0, 38513160,
            0, -2017356, 0, 92752, 0, -3927, 0, 210, 36, 2)
    )

    D = (
        1, 2, 6, 4, 30, 12, 42, 24, 90, 20, 66, 24, 2730, 420, 90, 48, 510,
        180, 3990, 840, 6930, 660, 690, 720, 13650, 1092, 378, 56, 870, 60,
        14322, 7392, 117810, 7140, 210, 72
    )

    pc = 36

    T = Type_IOP(p, pc, "pascal_triangle($p)"; msg)

    p ≥ 0 || throw(DomainError(p))
    
    o = p ≤ pc ? N[p] .// T(D[p]) : _faulhaber_BigInt(p)
    
    return o

end

# ------------------------------------------------------------------------------
#           faulhaber_polynomial(n::Integer, p::Int [; msg=true])
# ------------------------------------------------------------------------------

@doc raw"""
    faulhaber_polynomial(n::Integer, p::Int [; msg=true])

Faulhaber polynomial of degree `p` 
```math
    F(n,p)=\sum_{j=0}^{p}c_{j}n^{j},
```
where `n` is a positive integer and the coefficients are contained in the 
vector ``c=[c_0,⋯\ c_p]`` given by [`faulhaber_polynom`](@ref).

- `msg` : integer-overflow protection (IOP) - warning on activation
#### Examples:
```
julia> faulhaber_polynomial(3, 6)
276

julia> faulhaber_polynomial(5, 30)
IOP capture: faulhaber_polynomial(5, 30) autoconverted to Rational{BigInt}
186552813930161650665
```
"""
function faulhaber_polynomial(n::Integer, p::Int; msg=true)

    if n ≤ 0
        n < 0 ? throw(DomainError(n)) : return typeof(n)(0)
    elseif p ≤ 0
        p < 0 ? throw(DomainError(p)) : return typeof(n)(0)
    elseif p ≤ 36
        T = (float(n)^p < 9.223372036854776e12) ? typeof(n) : BigInt
    else
        T = BigInt
        str = "IOP - capture: "
        str *= "faulhaber_polynomial($n, $p) converted to Rational{BigInt}"
        msg && typeof(n) == Int && println(str)
    end

    F = CamiMath.faulhaber_polynom(T(p); msg=false)
    o = CamiMath.polynomial(F, T(n))

    Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")

    return Base.numerator(o)

end

# ------------------------------------------------------------------------------
#               faulhaber_summation(n::Integer, p::Int [; msg=true])
# ------------------------------------------------------------------------------

@doc raw"""
    faulhaber_summation(n::Integer, p::Int [; msg=true])

Sum of the ``p^{th}`` power of the first ``n`` natural numbers
```math
    \sum_{k=1}^{n}k^{p}=H_{n,-p}=F(n,p+1).
```
where ``H_{n,-p}`` is a [`harmonicNumber`](@ref)  of power `-p` and ``F(n,p)`` 
a [`faulhaber_polynomial`](@ref) of power `p`.

- `msg` : integer-overflow protection (IOP) - warning on activation
#### Examples:
```
julia> faulhaber_summation(3,5)
276

julia> faulhaber_summation(3,60)
IOP capture: faulhaber_polynom autoconverted to Rational{BigInt}
42391158276369125018901280178
```
"""
function faulhaber_summation(n::Integer, p::Int; msg=true)

    o = CamiMath.faulhaber_polynomial(n, p + 1; msg)

    return o

end