# SPDX-License-Identifier: MIT

# ==============================================================================
#                            Exponetial.jl
#                      Jook Walraven - 7-2-2023
# ==============================================================================

# ------------------------------------------------------------------------------
#               texp(x::T, a::T, p::Int) where {T<:Real}
# ------------------------------------------------------------------------------

@doc raw"""
    texp(x::T, a::T, p::Int) where T <: Real

Truncated exponential: Taylor expansion of ``exp(x)`` about ``x = a`` 
up to order `p``,
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
function texp(x::T, a::T, p::Int) where {T<:Real}

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

# ------------------------------------------------------------------------------
#                     log10_characteristic_power(x)
# ------------------------------------------------------------------------------

"""
    log10_characteristic_power(x)

characteristic power-of-10 of the number x
#### Examples:
```
log10_characteristic_power.([3,30,300])
3-element Vector{Int64}:
 0
 1
 2
```
"""
log10_characteristic_power(x) = Base.round(Int, Base.floor(log10(x)))

# ------------------------------------------------------------------------------
#                     log10_mantissa(x)
# ------------------------------------------------------------------------------

"""
    log10_mantissa(x)

log10 mantissa of the number x
#### Examples:
```
log10_mantissa.([3,30,300])
3-element Vector{Float64}:
 0.47712125471966244
 0.4771212547196624
 0.4771212547196626
```
"""
log10_mantissa(x) = Base.log10(x)-Base.floor(Base.log10(x))