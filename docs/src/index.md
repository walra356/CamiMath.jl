# CamMath.jl

Mathematics library with integer-overload protection

---
## Table of contents

```@contents
```

## Bernoulli number

```@docs
bernoulliB(n::T; msg=true) where {T<:Integer}   
bernoulliB_array(nmax::T; msg=true) where {T<:Integer}
```
## Factorial

```@docs
bigfactorial(n::T; msg=true) where {T<:Integer}
```

## Faulhaber polynomial

```@docs
faulhaber_polynom(k::Int; T=Int)
faulhaber_summation(n::Int, p::Int)
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