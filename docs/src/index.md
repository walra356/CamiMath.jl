# CamiMath.jl

Mathematics library with integer-overload protection

---
## Table of contents

```@contents
```

## Bernoulli number

```@docs
bernoulliB(n::Integer; msg=true, arr=false)
```
## Factorial

```@docs
bigfactorial(n::T; msg=true) where {T<:Integer}
```

## Faulhaber polynomial

```@docs
faulhaber_polynom(p::Integer; msg=true)
faulhaber_polynomial(n::Integer, p::Int; msg=true)
faulhaber_summation(n::Integer, p::Int)
```
## HarmonicNumber

```@docs
harmonicNumber(n::Integer; msg=true)
harmonicNumber(n::Integer, p::Int; msg=true)
```

## Pascal triangle

```@docs
pascal_triangle(row::Integer; msg=true)
pascal_next_row(a::Vector{Integer})
```