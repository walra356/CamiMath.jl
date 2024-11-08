# SPDX-License-Identifier: MIT

# ==============================================================================
#                            Polynomials.jl
#                    Jook Walraven - 19-2-2023
# ==============================================================================

# ------------------------------------------------------------------------------
#          polynomial(coords::NTuple{}, x::T; deriv=0) where {T<:Real}
# ------------------------------------------------------------------------------

function _polynom_primitive(coords)

    d = [1 // p for p ∈ Base.eachindex(coords)]

    coords = coords .* d

    return Base.pushfirst!(coords, 0)    # constant of integration equal to zero

end

# ..............................................................................
function _polynom_derivative(coords)

    d = Base.length(coords) - 1

    return [coords[p+1] * p for p = 1:d]

end

# ..............................................................................
function _polynom_derivatives(coords; deriv=0)

    for k = 1:deriv
        coords = _polynom_derivative(coords)
    end

    return coords

end

# ..............................................................................
@doc raw"""
    polynomial(coords, x::T [; deriv=0]) where T<:Real

Polynomial of degree ``d``,
```math
    P(x)=c_0 + c_1 x + ⋯ + c_d x^d,
```
where the coefficients `coords` = ``(c_0,⋯\ c_d)`` define the vector 
representation of the polynomial in a vector space of dimension ``d+1``.
#### Examples:
```
julia> coords = (1, 1, 1, 1, 1);
           
julia> P(x) = polynomial(coords,x);

julia> P(1)
5

julia> polynomial(coords, 1; deriv=1)     # P′(1)
10

julia> polynomial(coords, 2; deriv=2)     # P″(1)
20

julia> polynomial(coords,x; deriv=-1)   # primitive (zero integration constant)
137 // 60
```
"""
function polynomial(coords, x::T; deriv=0) where T<:Real

    isinteger(deriv) || error("Error: deriv not integer")

    d = Base.length(coords) - 1               # degree of polynomial

    coords = deriv == 0 ? coords :
             deriv ≥ d + 1 ? T(0) :
             deriv > 0 ? _polynom_derivatives(coords; deriv) :
             deriv == -1 ? _polynom_primitive(coords) :
             throw(DomainError(deriv))

    x === 1 && return sum(coords)

    d -= deriv
    o = coords[1]
    a = T(1)

    for i = 1:d
        a *= x
        o += (coords[1+i] * a)
    end

    return o

end

# ------------------------------------------------------------------------------
#                 polynom_power(coords, p::Int)
# ------------------------------------------------------------------------------

@doc raw"""
    polynom_power(coords, p)

The `coords` of a polynomial of degree ``d`` raised to the power `p`,
which define a polynomial in a vector space of dimension ``p d + 1``.
#### Examples:
```
julia> coords = (1,1,1)    # coordinates of polynomial vector of degree d = 2
(1, 1, 1)

julia> coords = (1,1,1);

julia> polynom_power(coords, 3)
7-element Vector{Int64}:
 1
 3
 6
 7
 6
 3
 1
```
"""
function polynom_power(coords, p::Int)

    p >= 0 || error("Error: negative polynom powers not allowed")
    p == 2 && return CamiMath.polynom_product(coords, coords)
    p == 1 && return coords
    p == 0 && return [1]

    o = CamiMath.polynom_product(coords, coords)

    for i = 1:p-2
        o = CamiMath.polynom_product(o, coords)
    end

    return o

end

# ------------------------------------------------------------------------------
#                    polynom_product(coords1, coords2)
# ------------------------------------------------------------------------------

@doc raw"""
    polynom_product(coords1, coords2)

The coefficients of the product of two polynomials, a = `coords1` and 
b = `coords2` of degree ``m`` and ``n``, which represents a polynomial 
in a vector space of dimension ``d=m+n+1``,
```math
    P(x)=a_0b_0 + (a_0b_1 + b_0a_1)x + ⋯ + a_n b_m x^{n+m}.
```
#### Examples:
```
julia> polynom_product((1, 1), (1, -1, 2))
4-element Vector{Int64}:
 1
 0
 1
 2

julia> polynom_product((1, 1), (1, -1.0, 2))
4-element Vector{Real}:
 1
 0.0
 1.0
 2  
```
"""
function polynom_product(coords1, coords2)

    a = coords1
    b = coords2

    n = Base.length(a)
    m = Base.length(b)

    if m > n

        o = [sum(a[1+j-i] * b[i] for i = 1:j) for j = 1:n]

        append!(o, [sum(a[n-i] * b[1+i+j] for i = 0:n-1) for j = 1:m-n])
        append!(o, [sum(a[n-i] * b[1+i+j+m-n] for i = 0:n-1-j) for j = 1:n-1])

    elseif m == n

        o = [sum(a[1+j-i] * b[i] for i = 1:j) for j = 1:n]

        append!(o, [sum(a[n-i] * b[1+i+j+m-n] for i = 0:n-1-j) for j = 1:n-1])

    else

        o = [sum(b[1+j-i] * a[i] for i = 1:j) for j = 1:m]

        append!(o, [sum(b[m-i] * a[1+i+j] for i = 0:m-1) for j = 1:n-m])
        append!(o, [sum(b[m-i] * a[1+i+j+n-m] for i = 0:m-1-j) for j = 1:m-1])

    end

    return o

end

# ------------------------------------------------------------------------------
#            polynom_product_expansion(coords1, coords2, p::Int)
# ------------------------------------------------------------------------------

@doc raw"""
    polynom_product_expansion(coords1, coords2, p::Int)

The coefficients of the product of two polynomials, a = `coords1` 
(of degree ``n``) and b = `coords2` (of degree ``m``), truncated at the 
order `p`, which represents a polynomial in a vector space of 
dimension ``d=p+1``
#### Examples:
```
julia> a = (1,-1,1);

julia> b = (1,1,-1,1,1,1);

julia> o = polynom_product(a, b); println(o)
[1, 0, -1, 3, -1, 1, 0, 1]
 
julia> o = polynom_product_expansion(a, b, 4); println(o)
[1, 0, -1, 3, -1] 
```
"""
function polynom_product_expansion(coords1, coords2, p::Int)

    a = coords1
    b = coords2

    n = Base.length(a)
    m = Base.length(b)
    n′ = n - 1
    m′ = m - 1

    if m > n

        o = [Base.sum(a[1+j-i] * b[1+i] for i = 0:j) for j = 0:min(n′, p)]
        p + 1 == length(o) && return o
        u = [sum(a[n-i] * b[1+i+j] for i = 0:n′) for j = 1:min(m - n, p - n′)]
        Base.append!(o, u)
        p + 1 == length(o) && return o

    elseif m == n

        o = [Base.sum(a[1+j-i] * b[1+i] for i = 0:j) for j = 0:min(n′, p)]
        p + 1 == length(o) && return o
        u = [sum(a[n-i] * b[i+j+m-n′] for i = 0:n′-j) for j = 1:min(n′, p - m′)]
        Base.append!(o, u)
        p + 1 == length(o) && return o

    else

        o = [Base.sum(b[1+j-i] * a[1+i] for i = 0:j) for j = 0:min(m′, p)]
        p + 1 == length(o) && return o
        u = [sum(b[m-i] * a[1+i+j] for i = 0:m′) for j = 1:min(n - m, p - m′)]
        Base.append!(o, u)
        p + 1 == length(o) && return o

    end

    return o

end