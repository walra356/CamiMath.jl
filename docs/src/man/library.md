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
isregular(sense::Type)
Type_IOP(n::Integer, nc::Integer, str="")
log10_characteristic(x)
log10_mantissa(x)
texp(x::T, a::T, p::Int) where {T<:Real}
```

## Bernoulli number

```@docs
bernoulliB(n::Integer; arr=false, msg=true)
```

## Divisor

```@docs
normalize_rationals(v::Vector{Rational{T}}) where {T<:Integer}
divisor(v::Vector{Rational{T}}) where {T<:Integer}
numerators(v::Vector{Rational{T}}) where {T<:Integer}
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
laguerreL(n::Integer, x::T; deriv=0, msg=true) where {T<:Real}
generalized_laguerreL(n::Integer, α, x::T; deriv=0, msg=true) where {T<:Real}
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
pochhammer(x::T, p::Int) where {T<:Real}
```

## Polynomials

```@docs
polynomial(coords, x::T; deriv=0) where {T<:Real}
polynom_power(coords, power::Int)
polynom_product(coords1, coords2)
polynom_product_expansion(coords1, coords2, p::Int)
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