# Toolbox

#### Dictionaries

```@docs
dictSuperscript
dictSubscript
dictUndoSmall
```

#### Superscript tools

```@docs
sup(i::T) where T<:Real
undo_sup(str::String)
```

#### Subscript tools

```@docs
sub(i::T) where T<:Real
undo_sub(str::String)
```

#### Tiny character tools

```@docs
undo_small(str::String)
frac(i::Rational{Int})
```

#### Rational covert to String

```@docs
strRational(n::T) where T<:Union{Rational{}, Int, BigInt}
```

## logical tools

```@docs
fwd
isforward(sense::Type)
bwd
isbackward(sense::Type)
reg
isregular(sense::Type)
rev
isreversed(sense::Type)
```

## Mathematics tools

```@docs
Type_IOP(n::Integer, nc::Integer, str="")
log10_characteristic(x)
log10_mantissa(x)
texp(x::T, a::T, p::Int) where T<:Real
```

## Divisor

```@docs
normalize_rationals(v::Vector{Rational{T}}) where T<:Integer
divisor(v::Vector{Rational{T}}) where T<:Integer
numerators(v::Vector{Rational{T}}) where T<:Integer
```
## Conversion to Big types

```@docs
convertToBig(x::T) where T
```
