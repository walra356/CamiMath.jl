# CamiMath.jl

Mathematics library with integer-overload protection (IOP)

---
## Table of contents

```@contents
```

## Bernoulli number

```@docs
bernoulliB(n::Integer; arr=false, msg=true)
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