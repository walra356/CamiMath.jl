# SPDX-License-Identifier: MIT

# ==============================================================================
#                          Laguerre.jl
#                    Jook Walraven - 1-3-2023
# ==============================================================================

# ------------------------------------------------------------------------------
#                  laguerre_polynom(n::Integer [; msg=true])
# ------------------------------------------------------------------------------

function _laguerre_polynom_BigInt(n, k)

    sgn = iseven(k) ? big(1) : -big(1)
    den = Base.factorial(big(n - k)) * Base.factorial(big(k))
    num = prod(k+1:n; init=sgn)

    return num .// den

end

@doc raw"""
    laguerre_polynom(n::Integer [; msg=true])
    
The coefficients of [`laguerreL`](@ref) for degree `n`, 
```math
    v_n=[c_0, c_1, \cdots\ c_n],
```
where, with ``k=0,1,⋯,n`` , 
```math
    c_k = \frac{\Gamma(n+1)}{\Gamma(k+1)}\frac{(-1)^{k}}{(n-k)!}\frac{1}{k!}.
```
- `msg` : integer-overflow protection (IOP) - warning on activation 
#### Example:
```
julia> laguerre_polynom(7)
(1//1, -7//1, 21//2, -35//6, 35//24, -7//40, 7//720, -1//5040)
```
"""
function laguerre_polynom(n::Integer; msg=true)

    N = (
        (1), (1, -1), (2, -4, 1), (6, -18, 9, -1), (24, -96, 72, -16, 1),
        (120, -600, 600, -200, 25, -1), (720, -4320, 5400, -2400, 450, -36, 1),
        (5040, -35280, 52920, -29400, 7350, -882, 49, -1),
        (40320, -322560, 564480, -376320, 117600, -18816, 1568, -64, 1),
        (362880, -3265920, 6531840, -5080320, 1905120, -381024, 42336, -2592,
            81, -1),
        (3628800, -36288000, 81648000, -72576000, 31752000, -7620480, 1058400,
            -86400, 4050, -100, 1),
        (39916800, -439084800, 1097712000, -1097712000, 548856000, -153679680,
            25613280, -2613600, 163350, -6050, 121, -1),
        (479001600, -5748019200, 15807052800, -17563392000, 9879408000,
            -3161410560, 614718720, -75271680, 5880600, -290400, 8712, -144, 1),
        (6227020800, -80951270400, 242853811200, -296821324800, 185513328000,
            -66784798080, 14841066240, -2120152320, 198764280, -12269400,
            490776, -12168, 169, -1),
        (87178291200, -1220496076800, 3966612249600, -5288816332800,
            3636061228800, -1454424491520, 363606122880, -59364264960,
            6492966480, -480960480, 24048024, -794976, 16562, -196, 1),
        (1307674368000, -19615115520000, 68652904320000, -99165306240000,
            74373979680000, -32724551059200, 9090153072000, -1669619952000,
            208702494000, -18036018000, 1082161080, -44717400, 1242150, -22050,
            225, -1),
        (20922789888000, -334764638208000, 1255367393280000, -1952793722880000,
            1586644899840000, -761589551923200, 232707918643200,
            -47491411968000, 6678479808000, -659602944000, 46172206080,
            -2289530880, 79497600, -1881600, 28800, -256, 1),
        (355687428096000, -6046686277632000, 24186745110528000,
            -40311241850880000, 35272336619520000, -18341615042150400,
            6113871680716800, -1372501805875200, 214453407168000,
            -23828156352000, 1906252508160, -110279070720, 4594961280,
            -135945600, 2774400, -36992, 289, -1),
        (6402373705728000, -115242726703104000, 489781588488192000,
            -870722823979008000, 816302647480320000, -457129482588979200,
            165074535379353600, -40426416827596800, 6948290392243200,
            -857813628672000, 77203226580480, -5104345559040, 248127909120,
            -8809274880, 224726400, -3995136, 46818, -324, 1)
    )

    D = (
        1, 1, 2, 6, 24, 120, 720, 5040, 40320, 62880, 3628800, 39916800,
        479001600, 6227020800, 87178291200, 1307674368000, 20922789888000,
        355687428096000, 6402373705728000
    )

    n ≥ 0 || throw(DomainError(n))

    nc = 18  # NB. length(D) = nc+1 (zero based notation)

    T = Type_IOP(n, nc; nam="laguerre_polynom", msg)

    if n ≤ nc
        return N[n+1] .// T(D[n+1])
    else
        return [_laguerre_polynom_BigInt(n, k) for k = 0:n]
    end

end

# ------------------------------------------------------------------------------
#                   _generalized_laguerre_coeff(n, α, m)
# ------------------------------------------------------------------------------

function _generalized_laguerre_Int(n::Int, α::Int, m::Int)

    sgn = iseven(m) ? 1 : -1
    den = Base.factorial(n - m) * factorial(m)
    num = prod(α+m+1:α+n; init=sgn)

    return num .// den

end

function _generalized_laguerre_BigInt(n::Integer, α::Integer, m::Integer)

    n = big(n)
    α = big(α)
    m = big(m)

    sgn = iseven(m) ? big(1) : -big(1)
    den = Base.factorial(n - m) * factorial(m)
    num = prod(α+m+big(1):α+n; init=sgn)

    return num .// den

end

function _generalized_laguerre_Float64(n::Int, α, m::Int)

    F = typeof(α)

    num = iseven(m) ? F(1.0) : -F(1.0)
    den = factorial(n - m) * factorial(m)
    den = F(den)

    for i = 1:(n-m)
        num *= F(α + m + i)
    end

    return num ./ den

end

function _generalized_laguerre_BigFloat(n::Integer, α, m::Integer)

    F = BigFloat

    num = iseven(m) ? F(1.0) : -F(1.0)
    den = factorial(big(n - m)) * factorial(big(m))
    den = F(den)

    for i = 1:(n-m)
        num *= F(α + m + i)
    end

    return num ./ den

end

# ------------------------------------------------------------------------------
#            generalized_laguerre_polynom(n::Int [, α=0 [; msg=true]])
# ------------------------------------------------------------------------------

@doc raw"""
    generalized_laguerre_polynom(n::Int [, α=0 [; msg=true]])

The coefficients of [`generalized_laguerreL`](@ref) for degree `n` and
parameter `α`,
```math
    _k(n, α) = \frac{\Gamma(α+n+1)}{\Gamma(α+k+1)}
    \frac{(-1)^{k}}{(n-k)!}\frac{1}{k!}
```
- `msg` : integer-overflow protection (IOP) - warning on activation 
#### Example:
```
julia> o =  generalized_laguerre_polynom(6,3); println(o)
Rational{Int64}[84//1, -126//1, 63//1, -14//1, 3//2, -3//40, 1//720]

julia> o =  generalized_laguerre_polynom(6,3.0); println(o)
[84.0, -126.0, 63.0, -14.0, 1.5, -0.075, 0.001388888888888889]
```
"""
function generalized_laguerre_polynom(n::Integer, α=0; msg=true)

    α === 0 && return laguerre_polynom(n; msg)

    n ≥ 0 || throw(DomainError(n))

    if α isa Integer
        T = Type_IOP(abs(α) + n, 18, α; nam="generalized_laguerre_polynom", msg)
        return T ≠ BigInt ?
               [_generalized_laguerre_Int(n, α, k) for k = 0:n] :
               [_generalized_laguerre_BigInt(n, α, k) for k = 0:n]
    else
        T = Type_IOP(n, 18, α; nam="generalized_laguerre_polynom", msg)
        return T ≠ BigInt ?
               [_generalized_laguerre_Float64(n, α, k) for k = 0:n] :
               [_generalized_laguerre_BigFloat(n, α, k) for k = 0:n]
    end

end

# ------------------------------------------------------------------------------
#              laguerreL(n::Int, α::U, x::T; deriv=0) where {U<:Real, T<:Real}
# ------------------------------------------------------------------------------

@doc raw"""
    laguerreL(n::Int, x::T; deriv=0) where T<:Real

Laguerre polynomal of degree `n`,
```math
    L_{n}(x)
    = \frac{1}{n!}e^{x}\frac{d^{n}}{dx^{n}}(e^{-x}x^{n})
    = \sum_{k=0}^{n}(-1)^{k}\binom{n}{n-k}\frac{x^{k}}{k!}
    = \sum_{k=0}^{n}c_k(n)x^{k}
```
where ``c_k(n)`` is a Laguerre coefficient from [`laguerre_polynom`](@ref).
#### Example:
```
julia> coords = laguerre_polynom(8); println(coords)
(1//1, -8//1, 14//1, -28//3, 35//12, -7//15, 7//180, -1//630, 1//40320)

julia> laguerreL(8, 5)
18029//8064

julia> polynomial(coords, 5)
18029//8064

julia> (xmin, Δx, xmax) = (0, 0.1, 11);
julia> n = 8;
julia> L = [laguerreL(n, x) for x=xmin:Δx:xmax];
julia> f = Float64.(L);
plot_function(f, xmin, Δx, xmax; title="laguerre polynomial (of degree $n)")
```
The plot is made using `CairomMakie`.
NB.: `plot_function` is not included in the `CamiXon` package.
![Image](../assets/laguerreL8.png)
"""
function laguerreL(n::Integer, x::T; deriv=0, msg=true) where {T<:Real}

    coords = laguerre_polynom(n; msg)

    o = polynomial(coords, x; deriv)

    return o

end



# ------------------------------------------------------------------------------
#              generalized_laguerreL(n::Int, α::U, x::T; deriv=0) where {U<:Real, T<:Real}
# ------------------------------------------------------------------------------

@doc raw"""
    generalized_laguerreL(n::Int, α::U, x::T; deriv=0) where {U<:Real, T<:Real}

Generalized Laguerre polynomal of degree `n` for parameter `α`,
```math
    L_{n}^{α}(x)
    = \frac{1}{n!}e^{x}x^{-α}\frac{d^{n}}{dx^{n}}(e^{-x}x^{n+α})
    = \sum_{k=0}^{n}(-1)^{k}\binom{n+α}{n-k}\frac{x^{k}}{k!}
    = \sum_{k=0}^{n}c_k(n,α)x^{k}
```
where ``c_k(n,α)`` is a generalized Laguerre coefficient from
[`generalized_laguerre_polynom`](@ref).
#### Example:
```
julia> coords = generalized_laguerre_polynom(5, 3); println(coords)
Rational{Int64}[56//1, -70//1, 28//1, -14//3, 1//3, -1//120]

julia> polynomial(coords, 10.0)
-10.666666666667311

julia> generalized_laguerreL(5, 10.0, 3)
-10.666666666667311
```
"""
function generalized_laguerreL(n::Integer, x::T, α=0; deriv=0, msg=true) where {T<:Real}

    coords = generalized_laguerre_polynom(n, α; msg)

    o = polynomial(coords, x; deriv)

    return o

end
