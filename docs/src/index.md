# CamiMath.jl

Mathematics library with integer-overload protection

---
## Table of contents

```@contents
```

## Bernoulli number

```@docs
bernoulliB(n::Integer; msg=true)  
bernoulliB_array(nmax::Integer; msg=true)
```
## Factorial

```@docs
bigfactorial(n::T; msg=true) where {T<:Integer}
```

## Faulhaber polynomial

```@docs
faulhaber_polynom(p::T; msg=true) where {T<:Integer}
faulhaber_polynomial(n::T, p::Int; msg=true) where {T<:Integer}
faulhaber_summation(n::T, p::Int) where {T<:Integer}
```
## HarmonicNumber

```@docs
harmonicNumber(n::T; msg=true) where {T<:Integer}
harmonicNumber_array(nmax::T; msg=true) where {T<:Integer}
harmonicNumber(n::T, p::Int; msg=true) where {T<:Integer}
```

## Pascal triangle

```@docs
pascal_triangle(nmax::T) where {T<:Integer}
pascal_next(a::Vector{Int})
```