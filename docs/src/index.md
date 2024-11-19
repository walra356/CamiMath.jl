```@meta
CurrentModule = CamiMath
```

# CamiMath.jl

Mathematics library with integer-overflow protection (IOP)

---

## Install

The package is installed using the Julia package manager

```
julia> using Pkg; Pkg.add("CamiMath")

julia> using CamiMath
```

## Table of contents

```@contents
```

## Singletons

```@docs
fwd
bwd
reg
rev
```

## Julia Toolbox

```@docs
isforward(sense::Type)
isbackward(sense::Type)
isregular(sense::Type)
isreversed(sense::Type)
Type_IOP(n::Integer, nc::Integer, str="")
log10_characteristic(x)
log10_mantissa(x)
texp(x::T, a::T, p::Int) where T<:Real
```

## Bernoulli number

```@docs
bernoulliB(n::Integer; arr=false, msg=true)
```

## Divisor

```@docs
normalize_rationals(v::Vector{Rational{T}}) where T<:Integer
divisor(v::Vector{Rational{T}}) where T<:Integer
numerators(v::Vector{Rational{T}}) where T<:Integer
```

## Factorial

```@docs
bigfactorial(n::Integer; msg=true)
```

## Faulhaber polynomial

```@docs
faulhaber_polynomial(n::Integer, p::Int; msg=true)
faulhaber_polynom(p::Integer; msg=true)
faulhaber_summation(n::Integer, p::Int)
```

## HarmonicNumber

```@docs
harmonicNumber(n::Integer; arr=false, msg=true)
harmonicNumber(n::Integer, p::Int; arr=false, msg=true)
```

## Fibonacci number

```@docs
fibonacci(n::Integer; arr=false, msg=true)
```

## Integer partitioning

```@docs
canonical_partitions(n::Int, m=0; header=true, reverse=true)
integer_partitions(n::Int, m=0; transpose=false, count=false)
```

## Lagrange polynomial

```@docs
lagrange_polynom(f::Vector{T}, start::Int, stop::Int, sense=fwd) where T<:Real
```

## Laguerre polynomial

```@docs
laguerreL(n::Integer, x::T; deriv=0, msg=true) where T<:Real
generalized_laguerreL(n::Integer, α, x::T; deriv=0, msg=true) where T<:Real
laguerre_polynom(p::Integer; msg=true)
generalized_laguerre_polynom(n::Integer, α=0; msg=true)
```

## Pascal triangle

```@docs
pascal_triangle(n::Integer; arr=false, msg=true)
pascal_next(a)
```

## Permutations

```@docs
permutations_unique_count(p::Vector{Vector{Int}}, i::Int)
```

## Pochhammer product

```@docs
pochhammer(x::T, p::Int) where T<:Real
```

## Polynomials

Polynomials can be regarded as the elements of vector spaces. As a first example we consider the set of all real polynomials of degree ``d``
```math
f_{\alpha}(x)=\alpha_{0}+\alpha_{1}x+\cdots+\alpha_{n}x^{d}.
```
These are maps 
```math
f_{\alpha}:\mathbb{\mathbb{R\rightarrow\mathbb{R}}} 
```
which satisfy the group operation addition of polynomials because the sum of two polynomials of degree ``d`` is again a polynomial of degree ``n``
```math
(f_{\alpha}+f_{\beta})(x)\equiv(\alpha_{0}+\beta_{0})+(\alpha_{1}+\beta_{1})x+\cdots+(\alpha_{n}+\beta_{n})x^{n}=f_{\alpha}(x)+f_{\beta}(x),
```
and remains a polynomial of degree n under scalar multiplication,
```math
(\lambda f_{\alpha})(x)\equiv f_{\lambda\alpha}(x)=\lambda\alpha_{0}+\lambda\alpha_{1}x+\cdots+\lambda\alpha_{n}x^{n}=\lambda f_{\alpha}(x).
```
The zero element of the vector space is the polynomial ``f_{\alpha}(x)=0`` and the inverse element of the element ``f_{\alpha}(x)`` is the polynomial ``-f_{\alpha}(x)``. Also the associative and distributive properties are easily verified. Hence, the set of all real polynomials of order n is a vector space over the field ``\mathbb{R}`` and the polynomials 
``1,x,x^{2},\cdots x^{d}`` represent a basis. This vector space is denoted by ``\mathcal{P}{}_{v}``. The coefficients ``\alpha_{0},\cdots\alpha_{v}`` are the coordinates of the vector ``f_{\alpha}`` with respect to this basis. Note that the set of all polynomials of degree ``m\leq d`` forms a subspace of the vector space of all polynomials of degree ``v``; i.e., ``\mathcal{P}{}_{m}\subseteq\mathcal{P}{}_{v}``. 

@doc raw"""
polynom :: Union{Vector{T}, Ntuple{T}} where T<:Real

The polynomial 
```
p(x) = c_0 + c_1 x + c_2 x^2 + ⋯ + c_3 x^p  
```
is defined by the coefficient vector
```
polynom = [c_0,c_1,⋯,c_p]
```
which represents the coordinates of a point in the vector space of the 
polynomial. 
the elements .
#### Example:
julia> polynom = [1,1,1,1,1]
[1,1,1,1,1]

julia> polynom = (1,1,1,1,1)
(1,1,1,1,1)
"""
struct polynom <: {Vector{T}, Tuple{T}} where T<:Real 
end

```@docs
polynomial(polynom, x::T; deriv=0) where T<:Real
polynom_power(polynom, power::Int)
polynom_product(polynom1, polynom2)
polynom_product_expansion(polynom1, polynom2, p::Int)
```

## Vector coupling

### Triangle relations

```@docs
istriangle(a::Real, b::Real, c::Real)
triangle_coefficient(a::Real, b::Real, c::Real)
```
### Vector-coupling coefficients

```@docs
threeJsymbol(j1::Real, m1::Real, j2::Real, m2::Real, j3::Real, m3::Real; msg=false)
CGC(j1::Real, m1::Real, j2::Real, m2::Real, J::Real, M::Real; msg=false)
```

## Index

```@index
```