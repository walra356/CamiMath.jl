## Bernoulli number

```@docs
bernoulliB(n::Integer; arr=false, msg=true)
```

## Divisor (common denominator of Rationals}

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
faulhaber_polynom(p::Integer; msg=true)
faulhaber_polynomial(n::Integer, p::Int; msg=true)
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

## Triangle relations

```@docs
istriangle(a::Real, b::Real, c::Real)
triangle_coefficient(a::Real, b::Real, c::Real)
```

## Truncated exponentials

```@docs
texp(x::T, a::T, p::Int) where {T<:Real}
log10_characteristic_power(x)
log10_mantissa(x)
```
