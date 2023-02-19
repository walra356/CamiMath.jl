# SPDX-License-Identifier: MIT

# ==============================================================================
#                            Triangle.jl
#                      Jook Walraven - 9-2-2023
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
Valid triangle condition: ``\Delta > 0``
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

    Δ = triangle_coefficient(a, b, c)

    o = Δ > 0 ? true : false

    return o

end