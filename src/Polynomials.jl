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
    polynomial(coords::NTuple{}, x::T [; deriv=0]) where {T<:Real}

The function ``f(x)=\text{polynomial}(c,x)``, where ``c=[c_0,⋯\ c_d]`` is the
the vector representation of a polynomial of degree `d` 
as given by [`polynom`](@ref).
```math
    \text{polynomial}(c,x)=c_0 + c_1 x + ⋯ + c_d x^d.
```
### Examples:
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
function polynomial(coords::NTuple{}, x::T; deriv=0) where {T<:Real}
    
    isinteger(deriv) || error("Error: deriv not integer")

    d = Base.length(coords) - 1               # degree of polynomial

    coords = deriv == 0 ? coords :
             deriv ≥ d + 1 ? T(0) :
             deriv > 0 ? _polynom_derivatives(coords; deriv) :
             deriv == -1 ? _polynom_primitive(coords) :
             throw(DomainError(deriv))

    d -= deriv
    v = Base.ones(T, d + 1)

    for i ∈ 1:d
        v[i+1] = v[i] * x
    end

    return sum(v .* coords)

end

# =============================== polynom_derivative_all(coords) =========

@doc raw"""
    polynom_derivatives_all(coords::Vector{<:Number})
Vector representation of all nontrivial derivatives of the polynomial `coords`.
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
### Examples:
```
coords=[1,1,1,1,1]               # vector representation of a polynomial of degree d=4
polynom_derivatives_all(coords)      # `all' (nontrivial) derivatives of polynomial `coords`
5-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [2, 6, 12]
 [6, 24]
 [24]
```
"""
function polynom_derivatives_all(coords::Vector{T}) where {T<:Real}

    k = Base.length(coords)

    coords = CamiXon.polynom_derivative(coords)

    o = [coords]

    for i = 2:k-1
        coords = polynom_derivative(coords)
        Base.push!(o, coords)
    end

    return o

end
# ==================================== polynom_power(coords, p) ================

@doc raw"""
    polynom_power(coords, p)
Vector representation of the polynomial `coords` raised to the power `p` which
results in a polynomial in a vector space of dimension ``p d + 1``.
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
### Examples:
```
coords=[1,1,1]             # vector representation of polynomial of degree ``d=2``
polynom_power(coords,2)
5-element Vector{Int64}:
 1
 2
 3
 2
 1
```
"""
function polynom_power(coords::Vector{T}, power::Int) where {T<:Real}

    power >= 0 || error("jwError: negative powers not allowed")
    power == 2 && return polynom_product(coords, coords)
    power == 1 && return coords
    power == 0 && return [1]

    o = CamiXon.polynom_product(coords, coords)

    for i = 1:power-2
        o = CamiXon.polynom_product(o, coords)
    end

    return o

end

# ==================================== polynom_powers(coords, pmax) ============

@doc raw"""
    polynom_powers(coords::Vector{T}, pmax::Int) where T<:Real
The polynomial `coords` raised to the powers 1,...,pmax  which
results in a collection of polynomials in vector spaces of dimension ``d+1`` tot ``p d + 1``.
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``.
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
### Examples:
```
coords=[1,1,1]                   # vector representation of polynomial of degree d=2
polynom_powers(coords,3)
3-element Vector{Vector{Int64}}:
 [1, 1, 1]
 [1, 2, 3, 2, 1]
 [1, 3, 6, 7, 6, 3, 1]
```
"""
function polynom_powers(coords::Vector{T}, pmax::Int) where {T<:Real}

    pmax > 0 || error("jwError: minimum power included is unity")

    o = [coords]

    for i = 1:pmax-1
        Base.push!(o, CamiXon.polynom_product(o[end], coords))
    end

    return o

end




# ================= polynom_product(a, b) ======================================

@doc raw"""
    polynom_product(a::Vector{T}, b::Vector{V}) where {T<:Real, V<:Real}
Vector representation of the product of two polynomials, ``a`` and ``b`` which
is a polynomial in a vector space of dimension ``d=m+n``,
```math
    p(c,x)=a_0b_0 + (a_0b_1 + b_0a_1)x + ⋯ + a_n b_m x^{n+m}.
```
Polynomials of degree ``d`` are represented by a vector in a vector space of dimension ``d+1``
The polynomial `coords` is specified by the coordinates vector ``c=[c_0,⋯\ c_d]``
consisting of the polynomial coefficients.
####
```
[polynom_product1([1.0,1],[1,-1,2])]
 [1.0, 0.0, 1.0, 2.0]
[polynom_product1([1//1,1],[1,-1,2])]
 [1//1, 0//1, 1//1, 2//1]
[polynom_product([1,1],[1,- 1,2])]
 [1, 0, 1, 2]
[polynom_product([1,- 1,2],[1,1])]
 [1, 0, 1, 2]
```
"""
function polynom_product(a::Vector{T}, b::Vector{V}) where {T<:Real,V<:Real}

    n = Base.length(a)
    m = Base.length(b)

    a, b = Base.promote(a, b)

    if m ≥ n
        o = [Base.sum(a[1+j-i] * b[1+i] for i = 0:j) for j = 0:n-1]
        if m ≠ n
            Base.append!(o, [Base.sum(a[n-i] * b[1+i+j] for i = 0:n-1) for j = 1:m-n])
        end
        Base.append!(o, [Base.sum(a[n-i] * b[1+i+j+m-n] for i = 0:n-1-j) for j = 1:n-1])
    else
        o = [Base.sum(b[1+j-i] * a[1+i] for i = 0:j) for j = 0:m-1]
        if m ≠ n
            Base.append!(o, [Base.sum(b[m-i] * a[1+i+j] for i = 0:m-1) for j = 1:n-m])
        end
        Base.append!(o, [Base.sum(b[m-i] * a[1+i+j+n-m] for i = 0:m-1-j) for j = 1:m-1])
    end

    return o

end

# ==================================== polynom_product_expansion(a, b, p) ============================================================

@doc raw"""
    polynom_product_expansion(a::Vector{T}, b::Vector{T}, p::Int) where T<:Real
Vector representation of the product of two polynomials, ``a`` (of degree ``n``) and ``b`` (of degree ``m``), with ``m≤n``
truncated at the order ``p`` is a polynomial in a vector space of dimension ``d=p+1``. If ``ab`` is the `polynom_product`,
the `polynom_product_expansion` is ``ab[1:p+1]``
####
```
a = [1,-1,1]
b = [1,1,-1,1,1,1]
o = polynom_product(a, b); println(o)
 [1, 0, -1, 3, -1, 1, 0, 1]
o = expand_product(a, b, 4); println(o)
 [1, 0, -1, 3, -1]
```
"""
function polynom_product_expansion(a::Vector{T}, b::Vector{T}, p::Int) where {T<:Real}

    n = Base.length(a)
    m = Base.length(b)

    if m ≥ n
        o = [Base.sum(a[1+j-i] * b[1+i] for i = 0:j) for j = 0:min(n - 1, p)]
        p + 1 == length(o) && return o
        if m ≠ n
            Base.append!(o, [Base.sum(a[n-i] * b[1+i+j] for i = 0:n-1) for j = 1:min(m - n, p - n + 1)])
        end
        p + 1 == length(o) && return o
        Base.append!(o, [Base.sum(a[n-i] * b[1+i+j+m-n] for i = 0:n-1-j) for j = 1:min(n - 1, p - m + 1)])
        p + 1 == length(o) && return o
    else
        o = [Base.sum(b[1+j-i] * a[1+i] for i = 0:j) for j = 0:min(m - 1, p)]
        p + 1 == length(o) && return o
        Base.append!(o, [Base.sum(b[m-i] * a[1+i+j] for i = 0:m-1) for j = 1:min(n - m, p - m + 1)])
        p + 1 == length(o) && return o
        Base.append!(o, [Base.sum(b[m-i] * a[1+i+j+n-m] for i = 0:m-1-j) for j = 1:min(m - 1, p - n + 1)])
        p + 1 == length(o) && return o
    end

    return o

end