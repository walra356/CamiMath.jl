# SPDX-License-Identifier: MIT

# ==============================================================================
#                          Faulhaber.jl
# ==============================================================================

global gl_faulhaber_Int = [
    0 // 1, 0 // 1, 90219075042845 // 12, 0 // 1, -593617720173709 // 24,
    0 // 1, 292938603384170 // 9, 0 // 1, -137675625585113 // 6,
    0 // 1, 10065214632473 // 1, 0 // 1, -18061762013765 // 6,
    0 // 1, 652976108950 // 1, 0 // 1, -429641365993 // 4,
    0 // 1, 249435438725 // 18, 0 // 1, -5758730593 // 4,
    0 // 1, 123028150 // 1, 0 // 1, -52802765 // 6, 0 // 1, 534905 // 1,
    0 // 1, -168113 // 6, 0 // 1, 11594 // 9, 0 // 1, -1309 // 24,
    0 // 1, 35 // 12, 1 // 2, 1 // 36
]

# ..............................................................................
global gl_faulhaber_BigInt = convert.(Rational{BigInt}, gl_faulhaber_Int)

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

# ..............................................................................
function _faulhaber_BigInt(p::T) where {T<:Integer}

    nul = big(0)
    one = big(1)

    p = big(p)

    P = CamMath.pascal_triangle(p)[end][1:end-1]
    B = CamMath.bernoulliB_array(p - one)  # was bernoulliB_array(p-1; T)
    B[2] = -B[2]

    F = (B .* P) // p

    F = Base.append!(F, nul // one) # add polynomial constant: c_0=0

    return Base.reverse(F)     # reverse to standard order

end

function faulhaber_polynom1(p::T) where {T<:Integer}

    str = "Warning: faulhaber_polynom converted to Rational{BigInt}"

    n = Int(p)
    nc = 36

    if n ≤ nc
        o = T == Int ? gl_faulhaber_Int[1:1+n] : gl_faulhaber_BigInt[1:1+n]
    else
        o = _faulhaber_BigInt(p)
        msg && T == Int && println(str)
    end

    return o

end
