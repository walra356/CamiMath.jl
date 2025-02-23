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
#                          Hermite.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                  hermite_polynom(n::Integer [; msg=true])
# ------------------------------------------------------------------------------

function _hermite_polynom_BigInt(n, k)

    B = BigInt

    if iseven(n+k) 
        m = (n-k)÷2
        sgn = iseven(m) ? B(1) : B(-1)
        num = prod(m+1:n; init=sgn) * B(2)^k 
        den = Base.factorial(B(k))
        o = num ÷ den
    else
        o = B(0)
    end

    return o

end

@doc raw"""
    hermite_polynom(n::Integer [; msg=true])
    
The coefficients of [`hermiteH`](@ref) for degree `n`, 
```math
    v(n)=[c_0, c_1, \cdots\ c_n],
```
where
```math
    c_k(n) = (-1)^{m}\frac{n!}{k!}\frac{2^{k}}{m!}
```
with  ``m = (n-k)/2`` and ``k=0,1,⋯,n``.

- `msg` : integer-overflow protection (IOP) - warning on activation 
#### Example:
```
julia> hermite_polynom(7)
(0, -1680, 0, 3360, 0, -1344, 0, 128)
```
"""
function hermite_polynom(n::Integer; msg=true)

    N = (
        (1), (0,2), (-2,0,4), (0,-12,0,8), (12,0,-48,0,16), (0,120,0,-160,0,32), (-120,0,720,0,-480,0,64), (0,-1680,0,3360,0,-1344,0,128), 
        (1680,0,-13440,0,13440,0,-3584,0,256),(0,30240,0,-80640,0,48384,0,-9216,0,512), (-30240,0,302400,0,-403200,0,161280,0,-23040,0,1024),
        (0,-665280,0,2217600,0,-1774080,0,506880,0,-56320,0,2048), (665280,0,-7983360,0,13305600,0,-7096320,0,1520640,0,-135168,0,4096),
        (0,17297280,0,-69189120,0,69189120,0,-26357760,0,4392960,0,-319488,0,8192), 
        (-17297280,0,242161920,0,-484323840,0,322882560,0,-92252160,0,12300288,0,-745472,0,16384),
        (0,-518918400,0,2421619200,0,-2905943040,0,1383782400,0,-307507200,0,33546240,0,-1720320,0,32768),
        (518918400,0,-8302694400,0,19372953600,0,-15498362880,0,5535129600,0,-984023040,0,89456640,0,-3932160,0,65536),
        (0,17643225600,0,-94097203200,0,131736084480,0,-75277762560,0,20910489600,0,-3041525760,0,233963520,0,-8912896,0,131072),
        (-17643225600,0,317578060800,0,-846874828800,0,790416506880,0,-338749931520,0,75277762560,0,-9124577280,0,601620480,0,-20054016,0,262144),
        (0,-670442572800,0,4022655436800,0,-6436248698880,0,4290832465920,0,-1430277488640,0,260050452480,0,-26671841280,0,1524105216,0,-44826624,0,524288)
    )

    n ≥ 0 || throw(DomainError(n))

    nc = 19 

    T = Type_IOP(n, nc; nam="hermite_polynom", msg)

    if n ≤ nc
        return N[n+1]
    else
        return [_hermite_polynom_BigInt(n, k) for k = 0:n]
    end

end

# ------------------------------------------------------------------------------
#       hermiteH(n::Integer, x::T; deriv=0, msg=true) where T<:Real
# ------------------------------------------------------------------------------

@doc raw"""
    hermiteH(n::Integer, x::T [; deriv=0 [, msg=true]]) where T<:Real

Hermite polynomal of degree `n`,
```math
    H_{n}(x)
    = (-1)^ne^{x^2}\frac{d^{n}}{dx^{n}}(e^{-x^2})
    = \sum_{k=0}^{n}(-1)^{m}\frac{n!}{k!}\frac{(2x)^{k}}{m!}
    = \sum_{k=0}^{n}c_k(n)x^{k}
```
where ``m = (n-k)/2`` and the ``c_k(n)`` are hermite coefficients of [`hermite_polynom`](@ref).
#### Example:
```
julia> hermite_polynom(7)
(0, -1680, 0, 3360, 0, -1344, 0, 128)

julia>  polynom = hermite_polynom(8)
(1680, 0, -13440, 0, 13440, 0, -3584, 0, 256)

julia> polynomial(polynom, 5)
52065680

julia> hermiteH(8, 5)
52065680
```
"""
function hermiteH(n::Integer, x::T; deriv=0, msg=true) where T<:Real

    polynom = hermite_polynom(n; msg)

    o = polynomial(polynom, x; deriv)

    return o

end