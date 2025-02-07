# Library

## Bernoulli number

```@docs
bernoulliB(n::Integer; arr=false, msg=true)
```

## Faulhaber polynomial

```@docs
faulhaber_polynomial(n::Integer, p::Int; msg=true)
faulhaber_polynom(p::Integer; msg=true)
faulhaber_summation(n::Integer, p::Int)
```

## Fibonacci number

```@docs
fibonacci(n::Integer; arr=false, msg=true)
```

## HarmonicNumber

```@docs
# harmonicNumber(n::Integer; arr=false, msg=true)
harmonicNumber(n::Integer, p::Int; arr=false, msg=true)
```

## Integer partitioning

```@docs
canonical_partitions(n::Int, m=0; header=true, reverse=true)
integer_partitions(n::Int, m=0; transpose=false, count=false)
```

## Lagrange polynom of tabulated function

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