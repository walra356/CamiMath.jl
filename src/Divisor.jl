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
#                            Divisor.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#        normalize_rationals(v::Vector{Rational{T}}) where {T<:Integer}
# ------------------------------------------------------------------------------

@doc raw"""
    normalize_rationals(v::Vector{Rational{T}}) where T<:Integer

Numerators separated from divisor
#### Example:
```
julia> normalize_rationals([1//1, 1//2, 1//3])
([6, 3, 2], 6)
```
"""
function normalize_rationals(v::Vector{Rational{T}}) where {T<:Integer}

    den = Base.gcd(v).den
    num = Base.convert(Vector{T}, (v .* den))

    return (num, den)

end

# ------------------------------------------------------------------------------
#        divisor(v::Vector{Rational{T}}) where {T<:Integer}
# ------------------------------------------------------------------------------

@doc raw"""
    divisor(v::Vector{Rational{T}}) where T<:Integer

Greatest common denominator of the set of rational numbers `v`
#### Example:
```
julia> divisor([1//1, 1//2, 1//3])
6
```
"""
function divisor(v::Vector{Rational{T}}) where T<:Integer

    return Base.gcd(v).den

end

# ------------------------------------------------------------------------------
#        numerators(v::Vector{Rational{T}}) where {T<:Integer}
# ------------------------------------------------------------------------------

@doc raw"""
    numerators(v::Vector{Rational{T}}) where T<:Integer

Numerators for the standard devisor of the set of rational numbers `v`
#### Example:
```
julia> numerators([1//1, 1//2, 1//3])
3-element Vector{Int64}:
 6
 3
 2
```
"""
function numerators(v::Vector{Rational{T}}) where T<:Integer

    den = Base.gcd(v).den
    num = Base.convert(Vector{T}, (v .* den))

    return num

end