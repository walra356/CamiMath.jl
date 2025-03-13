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
#                             Toolbox.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                            sup(i) 
# ------------------------------------------------------------------------------

function _superscript(i::Int)

    c = i < 0 ? [Char(0x207B)] : []

    for d ∈ reverse(digits(abs(i)))
        d == 0 ? Base.push!(c, Char(0x2070)) :
        d == 1 ? Base.push!(c, Char(0x00B9)) :
        d == 2 ? Base.push!(c, Char(0x00B2)) :
        d == 3 ? Base.push!(c, Char(0x00B3)) : Base.push!(c, Char(0x2070+d))
    end

    return join(c)

end
# ------------------------------------------------------------------------------
@doc raw"""
    sup(i::T) where T<:Real

Superscript notation for integers and rational numbers
#### Examples:
```
julia> sup(3) * 'P'
"³P"
```
"""
function sup(i::T) where T<:Real

    sgn = i < 0 ? Char(0x207B) : ""

    num = _superscript(numerator(abs(i)))
    den = _superscript(denominator(abs(i)))

    return T == Rational{Int} ? (sgn * num * '\U141F' * den) : (sgn * num)

end

# ------------------------------------------------------------------------------
#                                sub(i) 
# ------------------------------------------------------------------------------

function _subscript(i::Int)

    c = i < 0 ? [Char(0x208B)] : []

    for d ∈ reverse(digits(abs(i)))
        Base.push!(c, Char(0x2080+d))
    end

    return join(c)

end
# ------------------------------------------------------------------------------
function _subscript(str::String)

    d = Dict(
        'a' => Char(0x2090),
        'e' => Char(0x2091),
        'h' => Char(0x2095),
        'k' => Char(0x2096),
        'l' => Char(0x2097),
        'm' => Char(0x2098),
        'n' => Char(0x2099),
        'o' => Char(0x2092),
        'p' => Char(0x209A),
        'r' => Char(0x1D63),
        's' => Char(0x209B),
        't' => Char(0x209C),
        'x' => Char(0x2093))

    c::Vector{Char} = []

    for i ∈ collect(str)
        Base.push!(c, get(d, i,'.'))
    end

    return join(c)

end
# ------------------------------------------------------------------------------

@doc raw"""
    sub(i::T) where T<:Real
    sub(str::String)

Subscript notation for integers, rational numbers and a *subset* of lowercase 
characters ('a', 'e', 'h', 'k', 'l', 'm', 'n', 'o', 'p', 'r', 's', 't', 'x')
#### Examples:
```
julia> 'D' * sub(5//2)
"D₅⸝₂"

julia> "m" * sub("e")
"mₑ"
```
"""
function sub(i::T) where T<:Real

    sgn = i < 0 ? Char(0x208B) : ""

    num = _subscript(numerator(abs(i)))
    den = _subscript(denominator(abs(i)))

    return T == Rational{Int} ? (sgn * num * '\U2E1D' * den) : (sgn * num)

end
function sub(str::String)

    U = ['a','e','h','k','l','m','n','o','p','r','s','t','x']

    c = collect(str)

    for i ∈ eachindex(c)
        c[i] ∈ U || error("Error: subscript $(S[i]) not part of Unicode")
    end

    return _subscript(str::String)

end

# ------------------------------------------------------------------------------
#                                frac(i) 
# ------------------------------------------------------------------------------

@doc raw"""
    frac(i::Rational{Int})

Fraction notation for rational numbers
#### Examples:
```
julia> frac(-5//2)
"-⁵/₂"
```
"""
function frac(i::Rational{Int})

    sgn = i < 0 ? "-" : ""

    num = _superscript(numerator(abs(i)))
    den = _subscript(denominator(abs(i)))

    return sgn * num *  '/' * den

end

# ------------------------------------------------------------------------------
#                                strRational(i)
# ------------------------------------------------------------------------------

@doc raw"""
    strRational(n::T) where T<:Union{Rational{}, Int, BigInt}

Fraction notation for rational numbers and integers
#### Examples:
```
julia> strRational(-5//2)
"-5/2"
```
"""
function strRational(n::T) where T<:Union{Rational{}, Int, BigInt}


    isinteger(n) && return repr(n)

    sgn = n < 0 ? "-" : ""

    num = repr(numerator(abs(n)))
    den = repr(denominator(abs(n)))

    return sgn * num *  '/' * den

end

# ------------------------------------------------------------------------------
#                                 fwd
# ------------------------------------------------------------------------------
@doc raw"""
    fwd

Singleton type indicating `forward` sense
"""
struct fwd
end

# ------------------------------------------------------------------------------
#                            isforward(sense) 
# ------------------------------------------------------------------------------

@doc raw"""
    function isforward(sense)

Boolean status of `sense`, with options: [`fwd`](@ref) (forward) and [`bwd`](@ref) (backward).
#### Example:
```
julia> isforward(fwd)
true
```
"""
function isforward(sense)

strErr = "Error: invalid sense (options: fwd, bwd)"

return sense === fwd ? true : sense === bwd ? false : error(strErr)

end

# ------------------------------------------------------------------------------
#                                bwd
# ------------------------------------------------------------------------------
@doc raw"""
    bwd

Singleton type indicating `backward` sense
"""
struct bwd
end

# ------------------------------------------------------------------------------
#                            isbackward(sense) 
# ------------------------------------------------------------------------------

@doc raw"""
    function isbackward(sense)

Boolean status of `sense`, with options: [`fwd`](@ref) (forward) and [`bwd`](@ref) (backward).
#### Example:
```
julia> isbackward(fwd)
false
```
"""
function isbackward(sense)

strErr = "Error: invalid sense (options: fwd, bwd)"

return sense === bwd ? true : sense === fwd ? false : error(strErr)

end

# ------------------------------------------------------------------------------
#                                  reg
# ------------------------------------------------------------------------------
@doc raw"""
    reg

Singleton type indicating `regular` ordering
"""
struct reg
end

# ------------------------------------------------------------------------------
#                            isregular(sense) 
# ------------------------------------------------------------------------------

@doc raw"""
    function isregular(sense::Type)

Boolean status of `sense`, with options: [`reg`](@ref) (regular) and [`rev`](@ref) (reversed).
#### Example:
```
julia> isregular(reg)
true
```
"""
function isregular(sense)

strErr = "Error: invalid sense (options: reg, rev)" 

return sense === reg ? true : sense === rev ? false : error(strErr)

end

# ------------------------------------------------------------------------------
#                               rev
# ------------------------------------------------------------------------------

@doc raw"""
    rev

Singleton type indicating `reverse` ordering
"""
struct rev
end

# ------------------------------------------------------------------------------
#                            isreverse(sense) 
# ------------------------------------------------------------------------------

@doc raw"""
    function isreversed(sense::Type)

Boolean status of `sense`, with options: [`reg`](@ref) (regular) and [`rev`](@ref) (reversed).
#### Example:
```
julia> isreversed(rev)
true
```
"""
function isreversed(sense)

strErr = "Error: invalid sense (options: reg, rev)" 

return sense === rev ? true : sense === reg ? false : error(strErr)

end

# ============================= End ===========================

# ------------------------------------------------------------------------------
#                     log10_characteristic(x)   
# ------------------------------------------------------------------------------

"""
    log10_characteristic(x)

characteristic power-of-10 of the number `x`
#### Examples:
```
julia> log10_characteristic.([3,30,300])
3-element Vector{Int64}:
 0
 1
 2
```
"""
log10_characteristic(x) = Base.round(Int, Base.floor(log10(x)))

# ------------------------------------------------------------------------------
#                     log10_mantissa(x)
# ------------------------------------------------------------------------------

"""
    log10_mantissa(x)

log10 mantissa of the number `x`
#### Examples:
```
julia> log10_mantissa.([3,30,300])
3-element Vector{Float64}:
 0.47712125471966244
 0.4771212547196624
 0.4771212547196626
```
"""
log10_mantissa(x) = Base.log10(x) - Base.floor(Base.log10(x))

# ------------------------------------------------------------------------------
#         Type_IOP(n::Integer, nc::Integer [, a [; nam="" [; msg=true]]])
# ------------------------------------------------------------------------------

@doc raw"""
    Type_IOP(n::Integer, nc::Integer [, a [; nam="" [; msg=true]]])

`BigInt` if `n` is a `BigInt` or `n > nc`, otherwise `Int`; `a` is an 
auxiliary second variable.

- `nam` : function name

- `msg` : integer-overflow protection (IOP) - warning on activation 
#### Examples:
```
julia> Type_IOP(1, 1)
Int64

julia> Type_IOP(big(1), 1)
BigInt

julia> Type_IOP(2, 1)
BigInt

julia> Type_IOP(1, 1; nam="test")
Int64

julia> Type_IOP(2, 1, 0; nam="test")
 IOP capture at test(2, 0): output converted to BigInt
BigInt
```
"""
function Type_IOP(n::Integer, nc::Integer, a=nothing; nam="", msg=true)

    warning = " output converted to BigInt\n"

    if n isa BigInt
        return BigInt
    else
        isempty(nam) ? nothing :
        n ≤ nc ? nothing :
        !msg ? nothing :
        isnothing(a) ? print(" IOP capture at " * nam * "($n):" * warning) :
        print(" IOP capture at " * nam * "($n, $a):" * warning)
        return n ≤ nc ? Int : BigInt
    end

end

# ------------------------------------------------------------------------------
#               texp(x::T, a::T, p::Int) where T<:Real
# ------------------------------------------------------------------------------

@doc raw"""
    texp(x::T, a::T, p::Int) where T <: Real

Truncated exponential: Taylor expansion of ``exp(x)`` about ``x = a`` 
up to order `p`,
```math
    \mathsf{texp}(x,a,p) = 1+(x-a)+\frac{1}{2}(x-a)^2+⋯+\frac{1}{p!}(x-a)^p
```
### Examples:
```
julia> texp(1.0, 0.0, 5)
2.7166666666666663

julia> texp(1, 0, 5)
163//60
```
"""
function texp(x::T, a::T, p::Int) where T<:Real

    x = x - a

    o = y = T(1)

    x ≠ T(0) || return o

    if (T <: Rational) | (T <: Integer)
        for n = 1:p
            y *= x // T(n)
            o += y
        end
    else
        for n = 1:p
            y *= x / T(n)
            o += y
        end
    end

    return o

end

# ============================ convertToBig(x) ======================================

@doc raw"""
    convertToBig(x::T) where T

Conversion from 64-bit-based types to BigInt-based types:

* Int => BigInt,
* Vector{Int}                      => Vector{BigInt},
* Vector{Vector{Int}}              => Vector{Vector{BigInt}},
* Rational{Int}                    => Rational{BigInt},
* Vector{Rational{Int}}            => Vector{Rational{BigInt}},
* Vector{Vector{Rational{Int}}}    => Vector{Vector{Rational{BigInt}}},
* Float64                          => BigFloat,
* Vector{Float64}                  => Vector{BigFloat},
* Vector{Vector{Float64}}          => Vector{Vector{BigFloat}},
* Complex{Float64}                 => Complex{BigFloat},
* Vector{Complex{Float64}}         => Vector{Complex{BigFloat}},
* Vector{Vector{Complex{Float64}}} => Vector{Vector{Complex{BigFloat}}}

#### Example:
```
julia> [[1 // 1, 1 // 2], [1 // 1, 1 // 2]]
2-element Vector{Vector{Rational{Int64}}}:
 [1, 1//2]
 [1, 1//2]

julia> convertToBig([[1 // 1, 1 // 2], [1 // 1, 1 // 2]])
2-element Vector{Vector{Rational{BigInt}}}:
 [1, 1//2]
 [1, 1//2]
```
"""
function convertToBig(x::T) where T

    d = Dict(
        Int => BigInt,
        Vector{Int} => Vector{BigInt},
        Vector{Vector{Int}} => Vector{Vector{BigInt}},
        Rational{Int} => Rational{BigInt},
        Vector{Rational{Int}} => Vector{Rational{BigInt}},
        Vector{Vector{Rational{Int}}} => Vector{Vector{Rational{BigInt}}},
        Float64 => BigFloat,
        Vector{Float64} => Vector{BigFloat},
        Vector{Vector{Float64}} => Vector{Vector{BigFloat}},
        Complex{Float64} => Complex{BigFloat},
        Vector{Complex{Float64}} => Vector{Complex{BigFloat}},
        Vector{Vector{Complex{Float64}}} => Vector{Vector{Complex{BigFloat}}}
    )

    U = get(d, T, T)

    return convert(U, x)

end 