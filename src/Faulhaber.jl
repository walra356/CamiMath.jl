# SPDX-License-Identifier: MIT

# ==============================================================================
#                          Faulhaber.jl
# ==============================================================================

global gl_faulhaberInt = [
    [1 // 1], [0 // 1, 1 // 2, 1 // 2], [0 // 1, 1 // 6, 1 // 2, 1 // 3],
    [0 // 1, 0 // 1, 1 // 4, 1 // 2, 1 // 4], [0 // 1, -1 // 30, 0 // 1, 1 // 3,
        1 // 2, 1 // 5],
    [0 // 1, 0 // 1, -1 // 12, 0 // 1, 5 // 12, 1 // 2, 1 // 6],
    [0 // 1, 1 // 42, 0 // 1, -1 // 6, 0 // 1, 1 // 2, 1 // 2, 1 // 7],
    [0 // 1, 0 // 1, 1 // 12, 0 // 1, -7 // 24, 0 // 1, 7 // 12, 1 // 2,
        1 // 8],
    [0 // 1, -1 // 30, 0 // 1, 2 // 9, 0 // 1, -7 // 15, 0 // 1, 2 // 3, 1 // 2,
        1 // 9],
    [0 // 1, 0 // 1, -3 // 20, 0 // 1, 1 // 2, 0 // 1, -7 // 10, 0 // 1, 3 // 4,
        1 // 2, 1 // 10],
    [0 // 1, 5 // 66, 0 // 1, -1 // 2, 0 // 1, 1 // 1, 0 // 1, -1 // 1, 0 // 1,
        5 // 6, 1 // 2, 1 // 11],
    [0 // 1, 0 // 1, 5 // 12, 0 // 1, -11 // 8, 0 // 1, 11 // 6, 0 // 1,
        -11 // 8, 0 // 1, 11 // 12, 1 // 2, 1 // 12],
    [0 // 1, -691 // 2730, 0 // 1, 5 // 3, 0 // 1, -33 // 10, 0 // 1, 22 // 7,
        0 // 1, -11 // 6, 0 // 1, 1 // 1, 1 // 2, 1 // 13],
    [0 // 1, 0 // 1, -691 // 420, 0 // 1, 65 // 12, 0 // 1, -143 // 20, 0 // 1,
        143 // 28, 0 // 1, -143 // 60, 0 // 1, 13 // 12, 1 // 2, 1 // 14],
    [0 // 1, 7 // 6, 0 // 1, -691 // 90, 0 // 1, 91 // 6, 0 // 1, -143 // 10,
        0 // 1, 143 // 18, 0 // 1, -91 // 30, 0 // 1, 7 // 6, 1 // 2, 1 // 15],
    [0 // 1, 0 // 1, 35 // 4, 0 // 1,
        -691 // 24, 0 // 1, 455 // 12, 0 // 1, -429 // 16, 0 // 1, 143 // 12, 0 // 1, -91 // 24, 0 // 1,
        5 // 4, 1 // 2, 1 // 16], [0 // 1, -3617 // 510, 0 // 1, 140 // 3, 0 // 1, -1382 // 15, 0 // 1,
        260 // 3, 0 // 1, -143 // 3, 0 // 1, 52 // 3, 0 // 1, -14 // 3, 0 // 1, 4 // 3, 1 // 2, 1 // 17],
    [0 // 1, 0 // 1, -3617 // 60, 0 // 1, 595 // 3, 0 // 1, -11747 // 45, 0 // 1, 1105 // 6,
        0 // 1, -2431 // 30, 0 // 1, 221 // 9, 0 // 1, -17 // 3, 0 // 1, 17 // 12, 1 // 2, 1 // 18],
    [0 // 1, 43867 // 798, 0 // 1, -3617 // 10, 0 // 1, 714 // 1, 0 // 1, -23494 // 35, 0 // 1,
        1105 // 3, 0 // 1, -663 // 5, 0 // 1, 34 // 1, 0 // 1, -34 // 5, 0 // 1, 3 // 2, 1 // 2, 1 // 19],
    [0 // 1, 0 // 1, 43867 // 84, 0 // 1, -68723 // 40, 0 // 1, 2261 // 1, 0 // 1, -223193 // 140,
        0 // 1, 4199 // 6, 0 // 1, -4199 // 20, 0 // 1, 323 // 7, 0 // 1, -323 // 40, 0 // 1, 19 // 12,
        1 // 2, 1 // 20], [0 // 1, -174611 // 330, 0 // 1, 219335 // 63, 0 // 1, -68723 // 10, 0 // 1,
        6460 // 1, 0 // 1, -223193 // 63, 0 // 1, 41990 // 33, 0 // 1, -323 // 1, 0 // 1, 1292 // 21,
        0 // 1, -19 // 2, 0 // 1, 5 // 3, 1 // 2, 1 // 21], [0 // 1, 0 // 1, -1222277 // 220, 0 // 1,
        219335 // 12, 0 // 1, -481061 // 20, 0 // 1, 33915 // 2, 0 // 1, -223193 // 30, 0 // 1,
        146965 // 66, 0 // 1, -969 // 2, 0 // 1, 323 // 4, 0 // 1, -133 // 12, 0 // 1, 7 // 4, 1 // 2,
        1 // 22], [0 // 1, 854513 // 138, 0 // 1, -1222277 // 30, 0 // 1, 482537 // 6, 0 // 1,
        -755953 // 10, 0 // 1, 124355 // 3, 0 // 1, -223193 // 15, 0 // 1, 11305 // 3, 0 // 1,
        -3553 // 5, 0 // 1, 209 // 2, 0 // 1, -77 // 6, 0 // 1, 11 // 6, 1 // 2, 1 // 23],
    [0 // 1, 0 // 1, 854513 // 12, 0 // 1, -28112371 // 120, 0 // 1, 11098351 // 36, 0 // 1,
        -17386919 // 80, 0 // 1, 572033 // 6, 0 // 1, -5133439 // 180, 0 // 1, 37145 // 6, 0 // 1,
        -81719 // 80, 0 // 1, 4807 // 36, 0 // 1, -1771 // 120, 0 // 1, 23 // 12, 1 // 2, 1 // 24],
    [0 // 1, -236364091 // 2730, 0 // 1, 1709026 // 3, 0 // 1, -28112371 // 25, 0 // 1,
        22196702 // 21, 0 // 1, -17386919 // 30, 0 // 1, 208012 // 1, 0 // 1, -10266878 // 195,
        0 // 1, 29716 // 3, 0 // 1, -14421 // 10, 0 // 1, 506 // 3, 0 // 1, -253 // 15, 0 // 1, 2 // 1,
        1 // 2, 1 // 25], [0 // 1, 0 // 1, -1181820455 // 1092, 0 // 1, 21362825 // 6, 0 // 1,
        -28112371 // 6, 0 // 1, 277458775 // 84, 0 // 1, -17386919 // 12, 0 // 1, 1300075 // 3,
        0 // 1, -25667195 // 273, 0 // 1, 185725 // 12, 0 // 1, -24035 // 12, 0 // 1, 1265 // 6,
        0 // 1, -115 // 6, 0 // 1, 25 // 12, 1 // 2, 1 // 26], [0 // 1, 8553103 // 6, 0 // 1,
        -1181820455 // 126, 0 // 1, 55543345 // 3, 0 // 1, -52208689 // 3, 0 // 1,
        3606964075 // 378, 0 // 1, -20548177 // 6, 0 // 1, 2600150 // 3, 0 // 1, -10266878 // 63,
        0 // 1, 142025 // 6, 0 // 1, -16445 // 6, 0 // 1, 16445 // 63, 0 // 1, -65 // 3, 0 // 1,
        13 // 6, 1 // 2, 1 // 27], [0 // 1, 0 // 1, 76977927 // 4, 0 // 1, -3545461365 // 56, 0 // 1,
        166630035 // 2, 0 // 1, -469878201 // 8, 0 // 1, 721392815 // 28, 0 // 1, -61644531 // 8,
        0 // 1, 1671525 // 1, 0 // 1, -15400317 // 56, 0 // 1, 142025 // 4, 0 // 1, -29601 // 8,
        0 // 1, 4485 // 14, 0 // 1, -195 // 8, 0 // 1, 9 // 4, 1 // 2, 1 // 28], [0 // 1,
        -23749461029 // 870, 0 // 1, 179615163 // 1, 0 // 1, -709092273 // 2, 0 // 1,
        333260070 // 1, 0 // 1, -365460823 // 2, 0 // 1, 65581165 // 1, 0 // 1, -33193209 // 2,
        0 // 1, 3120180 // 1, 0 // 1, -905901 // 2, 0 // 1, 52325 // 1, 0 // 1, -9867 // 2, 0 // 1,
        390 // 1, 0 // 1, -273 // 10, 0 // 1, 7 // 3, 1 // 2, 1 // 29], [0 // 1, 0 // 1,
        -23749461029 // 60, 0 // 1, 5208839727 // 4, 0 // 1, -6854558639 // 4, 0 // 1,
        4832271015 // 4, 0 // 1, -10598363867 // 20, 0 // 1, 1901853785 // 12, 0 // 1,
        -137514723 // 4, 0 // 1, 22621305 // 4, 0 // 1, -8757043 // 12, 0 // 1, 303485 // 4,
        0 // 1, -26013 // 4, 0 // 1, 1885 // 4, 0 // 1, -609 // 20, 0 // 1, 29 // 12, 1 // 2, 1 // 30],
    [0 // 1, 8615841276005 // 14322, 0 // 1, -23749461029 // 6, 0 // 1, 15626519181 // 2,
        0 // 1, -102818379585 // 14, 0 // 1, 8053785025 // 2, 0 // 1, -31795091601 // 22, 0 // 1,
        731482225 // 2, 0 // 1, -137514723 // 2, 0 // 1, 19959975 // 2, 0 // 1, -2304485 // 2,
        0 // 1, 216775 // 2, 0 // 1, -16965 // 2, 0 // 1, 1131 // 2, 0 // 1, -203 // 6, 0 // 1,
        5 // 2, 1 // 2, 1 // 31], [0 // 1, 0 // 1, 8615841276005 // 924, 0 // 1,
        -736233291899 // 24, 0 // 1, 161474031537 // 4, 0 // 1, -3187369767135 // 112,
        0 // 1, 49933467155 // 4, 0 // 1, -328549279877 // 88, 0 // 1, 22675948975 // 28,
        0 // 1, -4262956413 // 32, 0 // 1, 68751025 // 4, 0 // 1, -14287807 // 8, 0 // 1,
        6720025 // 44, 0 // 1, -175305 // 16, 0 // 1, 2697 // 4, 0 // 1, -899 // 24, 0 // 1, 31 // 12,
        1 // 2, 1 // 32], [0 // 1, -7709321041217 // 510, 0 // 1, 68926730208040 // 693, 0 // 1,
        -2944933167596 // 15, 0 // 1, 184541750328 // 1, 0 // 1, -2124913178090 // 21, 0 // 1,
        36315248840 // 1, 0 // 1, -101092086116 // 11, 0 // 1, 36281518360 // 21, 0 // 1,
        -4262956413 // 17, 0 // 1, 28947800 // 1, 0 // 1, -57151228 // 21, 0 // 1, 2337400 // 11,
        0 // 1, -70122 // 5, 0 // 1, 7192 // 9, 0 // 1, -124 // 3, 0 // 1, 8 // 3, 1 // 2, 1 // 33],
    [0 // 1, 0 // 1, -84802531453387 // 340, 0 // 1, 17231682552010 // 21, 0 // 1,
        -16197132421778 // 15, 0 // 1, 761234720103 // 1, 0 // 1, -2337404495899 // 7, 0 // 1,
        99866934310 // 1, 0 // 1, -21662589882 // 1, 0 // 1, 49887087745 // 14, 0 // 1,
        -15630840181 // 34, 0 // 1, 47763870 // 1, 0 // 1, -28575614 // 7, 0 // 1, 292175 // 1,
        0 // 1, -89001 // 5, 0 // 1, 19778 // 21, 0 // 1, -682 // 15, 0 // 1, 11 // 4, 1 // 2, 1 // 34],
    [0 // 1, 2577687858367 // 6, 0 // 1, -84802531453387 // 30, 0 // 1,
        117175441353668 // 21, 0 // 1, -78671786048636 // 15, 0 // 1, 2875775609278 // 1,
        0 // 1, -7224704805506 // 7, 0 // 1, 261190443580 // 1, 0 // 1, -245509351996 // 5,
        0 // 1, 49887087745 // 7, 0 // 1, -822675799 // 1, 0 // 1, 77331980 // 1, 0 // 1,
        -42242212 // 7, 0 // 1, 397358 // 1, 0 // 1, -336226 // 15, 0 // 1, 23188 // 21, 0 // 1,
        -748 // 15, 0 // 1, 17 // 6, 1 // 2, 1 // 35], [0 // 1, 0 // 1, 90219075042845 // 12, 0 // 1,
        -593617720173709 // 24, 0 // 1, 292938603384170 // 9, 0 // 1, -137675625585113 // 6,
        0 // 1, 10065214632473 // 1, 0 // 1, -18061762013765 // 6, 0 // 1, 652976108950 // 1,
        0 // 1, -429641365993 // 4, 0 // 1, 249435438725 // 18, 0 // 1, -5758730593 // 4, 0 // 1,
        123028150 // 1, 0 // 1, -52802765 // 6, 0 // 1, 534905 // 1, 0 // 1, -168113 // 6, 0 // 1,
        11594 // 9, 0 // 1, -1309 // 24, 0 // 1, 35 // 12, 1 // 2, 1 // 36]
]

# ..............................................................................
global gl_faulhaberBigInt = convert.(Vector{Rational{BigInt}}, gl_faulhaberInt)

# ..............................................................................
function _faulhaber_BigInt(p::T) where {T<:Integer}

    nul = big(0)
    one = big(1)

    p = big(p)

    P = CamMath.pascal_triangle(p)[end][1:end-1]
    B = CamMath.bernoulliB_array(p - one)  # was bernoulliB_array(p-1; T)
    B[2] = -B[2]

    F = (B .* P) // p

    F = Base.append!(F, nul // one)  # add polynomial constant: c_0 = 0//1

    return Base.reverse(F)     # reverse to standard order

end

# ..............................................................................
@doc raw"""
    faulhaber_polynom(p::T; msg=true) where {T<:Integer}

Vector representation of the coefficients of the [`faulhaber_polynomial`](@ref) 
of degree `p` 
```math
   c=[c_0,⋯\ c_p]
```
where
```math
    c_0=0,\ \ c_j=\frac{1}{p}{\binom{p}{p-j}}B_{p-j},
```
with ``j∈\{ 1,⋯\ p\}``. The ``B_0,⋯\ B_{p-1}`` are Bernoulli numbers
(but with ``B_1=+\frac{1}{2}`` rather than ``-\frac{1}{2}``).
Integer-overflow protection: for `p > 36` the output is autoconverted to 
Rational{BigInt}. By default the capture message is activated: 
"Warning: faulhaber_polynom converted to Rational{BigInt}". 
### Example:
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
function faulhaber_polynom(p::T; msg=true) where {T<:Integer}

    str = "Warning: faulhaber_polynom converted to Rational{BigInt}"

    n = Int(p)
    nc = 36

    if n ≤ nc
        o = T == Int ? gl_faulhaberInt[n][1:end] : gl_faulhaberBigInt[n][1:end]
    else
        o = _faulhaber_BigInt(p)
        msg && T == Int && println(str)
    end

    return o

end

# ==================== faulhaber_polynomial(n,p) ===============================

@doc raw"""
    faulhaber_polynomial(n::T, p::Int; msg=true) where {T<:Integer}

Faulhaber polynomial of degree `p` 
```math
    F(n,p)=\sum_{j=0}^{p}c_{j}n^{j},
```
where the coefficients are contained in the coefficient vector 
[`faulhaber_polynom`](@ref).
### Example:
```
julia> faulhaber_polynomial(3, 6)
276
```
"""
function faulhaber_polynomial(n::T, p::Int; msg=true) where {T<:Integer}

    str = "Warning: faulhaber_polynomial converted to Rational{BigInt}"

    n ≠ 0 || return 0

    if n^p < 9e18
        W = T
    else
        W = BigInt
        msg && T == Int && println(str)
    end

    F = CamMath.faulhaber_polynom(W(p))
    o = 0
    for k = 1:p
        for i = 1:k
            F[1+k] *= n # avoid n^k in o = Base.sum([F[k+1]*n^k for k=1:p+1])
        end
        o += F[1+k]
    end

    Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")

    return Base.numerator(o)

end

# =================================== faulhaber_summation(n,p;T) ===============

@doc raw"""
    faulhaber_summation(n::T, p::Int) where {T<:Integer}

Sum of the ``p^{th}`` power of the first ``n`` natural numbers
```math
    \sum_{k=1}^{n}k^{p}=F(n,p+1).
```
where ``F(n,p)`` is the [`faulhamer_polynomial`](@ref) of degree `p`.
### Examples:
```
julia> faulhaber_summation(3,5)
276

julia> faulhaber_summation(3,60)
Warning: faulhaber_polynom converted to Rational{BigInt}
42391158276369125018901280178
```
"""
function faulhaber_summation(n::T, p::Int) where {T<:Integer}

    o = faulhaber_polynomial(n, p + 1)

    return o

end