# SPDX-License-Identifier: MIT

# ==============================================================================
#                          Faulhaber.jl
# ==============================================================================

global gl_faulhaber_Int = [
    [1 // 1], [0//1, 1//2, 1//2], [0//1, 1//6, 1//2, 1//3],
    [0//1, 0//1, 1//4, 1//2, 1//4], [0//1, -1//30, 0//1, 1//3, 1//2, 1//5],
    [0//1, 0//1, -1//12, 0//1, 5//12, 1//2, 1//6],
    [0//1, 1//42, 0//1, -1//6, 0//1, 1//2, 1//2, 1//7],
    [0//1, 0//1, 1//12, 0//1, -7//24, 0//1, 7//12, 1//2, 1//8],
    [0//1, -1//30, 0//1, 2//9, 0//1, -7//15, 0//1, 2//3, 1//2, 1//9],
    [0//1, 0//1, -3//20, 0//1, 1//2, 0//1, -7//10, 0//1, 3//4, 1//2, 1//10],
    [0//1, 5//66, 0//1, -1//2, 0//1, 1//1, 0//1, -1//1, 0//1, 5//6, 1//2, 1//11],
    [0//1, 0//1, 5//12, 0//1, -11//8, 0//1, 11//6, 0//1, -11//8, 0//1, 11//12, 1//2, 1//12],
    [0//1, -691//2730, 0//1, 5//3, 0//1, -33//10, 0//1, 22//7, 0//1, -11//6, 0//1, 1//1, 1//2, 1//13],
    [0//1, 0//1, -691//420, 0//1, 65//12, 0//1, -143//20, 0//1, 143//28, 0//1, -143//60, 0//1, 13//12, 1//2, 1//14],
    [0//1, 7//6, 0//1, -691//90, 0//1, 91//6, 0//1, -143//10, 0//1, 143//18, 0//1, -91//30, 0//1, 7//6, 1//2, 1//15],
    [0//1, 0//1, 35//4, 0//1, -691//24, 0//1, 455//12, 0//1, -429//16, 0//1, 143//12, 0//1, -91//24, 0//1, 5//4, 1//2, 1//16],
    [0//1, -3617//510, 0//1, 140//3, 0//1, -1382//15, 0//1, 260//3, 0//1, -143//3, 0//1, 52//3, 0//1, -14//3, 0//1, 4//3, 1//2, 1//17],
    [0//1, 0//1, -3617//60, 0//1, 595//3, 0//1, -11747//45, 0//1, 1105//6, 0//1, -2431//30, 0//1, 221//9, 0//1, -17//3, 0//1, 17//12, 1//2, 1//18],
    [0//1, 43867//798, 0//1, -3617//10, 0//1, 714//1, 0//1, -23494//35, 0//1, 1105//3, 0//1, -663//5, 0//1, 34//1, 0//1, -34//5, 0//1, 3//2, 1//2, 1//19]
]

# ..............................................................................
global gl_faulhaber_BigInt = convert.(Vector{Rational{BigInt}}, gl_faulhaber_Int)

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
    faulhaber_polynom(p::T) where {T<:Integer}

Vector representation of the coefficients of the [`faulhaber_polynomial`](@ref) 
of degree `p` 
```math
   c=[c_0,⋯\ c_p]
```
with vector elements
```math
    c_0=0,\ \ c_j=\frac{1}{p}{\binom{p}{p-j}}B_{p-j},
```
where ``j∈\{ 1,⋯\ p\}``. The ``B_0,⋯\ B_{p-1}`` are Bernoulli numbers
(but with ``B_1=+\frac{1}{2}`` rather than ``-\frac{1}{2}``).
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
function faulhaber_polynom1(p::T; msg=true) where {T<:Integer}

    str = "Warning: faulhaber_polynom converted to Rational{BigInt}"

    n = Int(p)
    nc = 36

    if n ≤ nc
        o = T == Int ? gl_faulhaber_Int[n][1:end] : gl_faulhaber_BigInt[n][1:end]
    else
        o = _faulhaber_BigInt(p)
        msg && T == Int && println(str)
    end

    return o

end

# ==================== faulhaber_polynomial(n,p) ===============================

@doc raw"""
    faulhaber_polynomial(n::T, p::T) where {T<:Integer}

Faulhaber polynomial of degree `p` 
```math
    F(n,p)=\sum_{j=0}^{p}c_{j}n^{j},
```
where the coefficients are contained in the coefficient vector 
[`faulhaber_polynom`](@ref).
### Example:
```
```
"""
function faulhaber_polynomial(n::T, p::T) where {T<:Integer}

    n ≠ 0 || return 0

    F = CamMath.faulhaber_polynom(p)
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
    faulhaber_summation(n::T, p::T) where {T<:Integer}

Sum of the ``p^{th}`` power of the first ``n`` natural numbers
```math
    \sum_{k=1}^{n}k^{p}=F(n,p+1).
```
where ``F(n,p)`` is the [`faulhamer_polynomial`](@ref) of degree `p`.
### Examples:
```
faulhaber_summation(5,1)
 15

faulhaber_summation(3,60; T=BigInt)
  42391158276369125018901280178
```
"""
function faulhaber_summation(n::T, p::T) where {T<:Integer}

    o = faulhaber_polynomial(n, p + 1)

    return o

end

function faulhaber_polynom(p::T) where {T<:Integer}

    p < 1 && return 0
    p > 1 || return 1 // 1

    P = CamMath.pascal_triangle(p)[end][1:end-1]
    B = CamMath.bernoulliB_array(p - 1)  # was bernoulliB_array(p-1; T)
    B[2] = -B[2]

    F = (B .* P) // p

    F = Base.append!(F, 0 // 1)   # add polynomial constant (zero in this case)

    return Base.reverse(F)     # reverse to standard order

end