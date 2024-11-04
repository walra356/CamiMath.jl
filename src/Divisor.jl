# SPDX-License-Identifier: MIT

# ==============================================================================
#                            Divisor.jl
#                      Jook Walraven - 18-2-2023
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