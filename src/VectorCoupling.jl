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
#                         AngularMomentum.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#            triangle_coefficient(a::Real, b::Real, c::Real)
# ------------------------------------------------------------------------------

@doc raw"""
    triangle_coefficient(a::Real, b::Real, c::Real)

Triangle coefficient for a triangle of sides `a`, `b` and `c`,
```math
    \Delta(abc)\equiv\frac{(a+b-c)!(b+c-a)!(c+a-b)!}{(a+b+c+1)!}
```
Triangle condition satisfied for ``\Delta > 0``
#### Examples:
```
julia> triangle_coefficient(3, 4, 5)
1//180180

julia> triangle_coefficient(1//2, 1, 1.5)
1//12
```
"""
function triangle_coefficient(a::Real, b::Real, c::Real)

    (a, b, c) = Base.promote(a, b, c)

    Base.isinteger(a + b + c) || return 0

    A = Int(a + b - c)
    B = Int(b + c - a)
    C = Int(c + a - b)

    A = A ≥ 0 ? CamiMath.bigfactorial(A) : return 0
    B = B ≥ 0 ? CamiMath.bigfactorial(B) : return 0
    C = C ≥ 0 ? CamiMath.bigfactorial(C) : return 0

    N = A * B * C
    D = CamiMath.bigfactorial(Int(a + b + c + 1))

    return N // D

end

# ------------------------------------------------------------------------------
#               istriangle(a::Real, b::Real, c::Real)
# ------------------------------------------------------------------------------

@doc raw"""
    istriangle(a::Real, b::Real, c::Real)

Triangle condition for a triangle of sides `a`, `b` and `c`.
#### Example:
```
julia> istriangle(3, 4, 5)
true

julia> istriangle(1//2, 1, 1.5)
true
```
"""
function istriangle(a::Real, b::Real, c::Real)

    Δ = CamiMath.triangle_coefficient(a, b, c)

    o = Δ > 0 ? true : false

    return o

end

# ------------------------------------------------------------------------------
#              threeJsymbol(j1, m1, j2, m2, j3, m3 [; msg=true])
# ------------------------------------------------------------------------------

function _strRational(n::T) where {T<:Union{Rational{},Int,BigInt}}


    isinteger(n) && return repr(n)

    sgn = n < 0 ? "-" : ""

    num = repr(numerator(abs(n)))
    den = repr(denominator(abs(n)))

    return sgn * num * '/' * den

end
# ........................................................
function _Racah_sqrt2(j1, m1, j2, m2, J, M)

    a = Int(j1 + m1)
    b = Int(j1 - m1)
    c = Int(j2 + m2)
    d = Int(j2 - m2)
    e = Int(J + M)
    f = Int(J - M)

    a = a ≥ 0 ? factorial(big(a)) : return 0
    b = b ≥ 0 ? factorial(big(b)) : return 0
    c = c ≥ 0 ? factorial(big(c)) : return 0
    d = d ≥ 0 ? factorial(big(d)) : return 0
    e = e ≥ 0 ? factorial(big(e)) : return 0
    f = f ≥ 0 ? factorial(big(f)) : return 0

    return a * b * c * d * e * f

end
# ........................................................
function _Racah_denom(j1, m1, j2, m2, J, t::Int)

    a = Int(J - j2 + t + m1)
    b = Int(J - j1 + t - m2)
    c = Int(j1 + j2 - J - t)
    d = Int(j1 - t - m1)
    e = Int(j2 - t + m2)

    a = a ≥ 0 ? factorial(big(a)) : return 0
    b = b ≥ 0 ? factorial(big(b)) : return 0
    c = c ≥ 0 ? factorial(big(c)) : return 0
    d = d ≥ 0 ? factorial(big(d)) : return 0
    e = e ≥ 0 ? factorial(big(e)) : return 0

    return a * b * c * d * e * factorial(big(t))

end
# ........................................................
function _Racah_sum(j1, m1, j2, m2, J)

    o = big(0)

    for t = 0:(j1+j2-J)
        sign = iseven(t) ? 1 : -1
        d = _Racah_denom(j1, m1, j2, m2, J, t)
        o += d > 0 ? sign // d : 0
    end

    return o

end
# ........................................................

@doc raw"""
    threeJsymbol(j1::Real, m1::Real, j2::Real, m2::Real, j3::Real, m3::Real [; msg=true])

Wigner 3j symbol. This is a vector coupling coefficient with optimized symmetry
properties. The 3j symbols are zero unless ``Δ(j_{1},j_{2},j_{3})>0``
(triangle inequality holds) and ``m_{1}+m_{2}+m_{3}=0``. The implementation
is based on the Racah formula:

```math
\left(\begin{array}{ccc}
j_{1} & j_{2} & j_{3}\\
m_{1} & m_{2} & m_{3}
\end{array}\right)=
(-1)^{j_{1}-j_{2}-m_{3}}\sqrt{\Delta(j_{1}j_{2}J)}\\\times
\sqrt{\left(j_{1}+m_{1}\right)!
\left(j_{1}-m_{1}\right)!
\left(j_{2}+m_{2}\right)!
\left(j_{2}-m_{2}\right)!
\left(j_{3}+m_{3}\right)!
\left(j_{3}-m_{3}\right)!}
\\\times\sum_{t}\frac{(-)^{t}}{t!(j_{3}-j_{2}+t+m_{1})!
(j_{3}-j_{1}+t-m_{2})!
(j_{1}+j_{2}-j_{3}-t)!(j_{1}-t-m_{1})!(j_{2}-t+m_{2})!}
```
#### Example:
```
julia> o = threeJsymbol(3, 0, 4, -1, 5, 1); println(" = $o")
-√(361/30030) = -0.1096417439724123565166029917781360897459044055433631161836138910409772907333476

julia> threeJsymbol(3, 0, 4, -1, 5, 1; msg=false)
-0.1096417439724123565166029917781360897459044055433631161836138910409772907333476

julia> threeJsymbol(0, 0, 0, 0, 0, 0; msg=false)
1.0
```
"""
function threeJsymbol(j1::Real, m1::Real, j2::Real, m2::Real, j3::Real, m3::Real; msg=true)

    (j1, m1, j2, m2, j3, m3) = promote(j1, m1, j2, m2, j3, m3)

    iszero(m1 + m2 + m3) || return 0

    Δ = CamiMath.triangle_coefficient(j1, j2, j3)
    T = _Racah_sqrt2(j1, m1, j2, m2, j3, m3)
    R = _Racah_sum(j1, m1, j2, m2, j3)
    S = R * R
    A = Δ * T * S

    sgn_phase = iseven(j1 - j2 + m3) ? 1 : -1
    sgn_racah = sign(R)
    sgn = sgn_phase * sgn_racah

    msg && print((sgn < 0 ? "-" : "") * "√(" * _strRational(A) * ")")

    return sgn * sqrt(A)

end

# ------------------------------------------------------------------------------
#              CGC(j1, m1, j2, m2l, J, M [; msg=true])
# ------------------------------------------------------------------------------

@doc raw"""
    CGC(j1::Real, m1::Real, j2::Real, m2::Real, J::Real, M::Real [; msg=true])

Clebsch-Gordan coefficient (CGC). This is a vector-coupling coefficient in
Dirac notation. The CGCs are zero unless ``Δ(j_{1},j_{2},j_{3})>0``
(triangle inequality holds) and ``M=m_{1}+m_{2}``. The relation to the
Wigner 3j symbols is given by:

```math
\langle j_{1}m_{1};j_{2}m_{2}|JM\rangle\equiv
(-1)^{j_{1}-j_{2}+M}\sqrt{2J+1}\left(\begin{array}{ccc}
j_{1} & j_{2} & J\\
m_{1} & m_{2} & -M
\end{array}\right)
```
#### Example:
```
julia> j1=3; m1=0; j2=4; m2=-1; J=5; M=-1;

julia> o = CGC(j1, m1, j2, m2, J, M); println(" = $o")
-√(361/2730) = -0.36364052611670255269921486774521555203216489725107181148303161368088211274565

julia> o = CGC(j1, m1, j2, m2, J, M; msg=false); println(o)
-0.36364052611670255269921486774521555203216489725107181148303161368088211274565

julia> o = (-1)^(j1-j2+M) * sqrt(2J+1) * threeJsymbol(j1, m1, j2, m2, J, -M; msg=false); println(o)
-0.36364052611670255269921486774521555203216489725107181148303161368088211274565

julia> CGC(3, 0, 4, -1, 5, -1);
-√(361/2730)
```
"""
function CGC(j1::Real, m1::Real, j2::Real, m2::Real, J::Real, M::Real; msg=true)

    (j1, m1, j2, m2, J, M) = promote(j1, m1, j2, m2, J, M)

    sgn = iseven(j1 - j2 + M) ? 1 : -1
    tJs = threeJsymbol(j1, m1, j2, m2, J, -M; msg=false)

    if msg
        Δ = CamiMath.triangle_coefficient(j1, j2, J)
        T = _Racah_sqrt2(j1, m1, j2, m2, J, M)
        R = _Racah_sum(j1, m1, j2, m2, J)
        S = R * R
        A = Δ * T * S
        s = sign(R) < 0 ? "-" : ""
        print(s * "√(" * _strRational(A * (2J + 1)) * ")")
    end

    return sgn * sqrt(2J + 1) * tJs

end
