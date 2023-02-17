# ====================== permutations_unique_count(p, i) =======================

@doc raw"""
    permutations_unique_count(p::Array{Array{Int64,1},1}, i::Int)
Number of unique permutations of the subarray ``p[i]``.
#### Example:
```
p = [[1,2,3],[2,3,1,4,3]]
permutations_unique_count(p,2)
 60
```
"""
function permutations_unique_count(p::Array{Array{Int64,1},1}, i::Int)

    o = Base.factorial(Base.length(p[i]))
    d = Base.Dict([(n,Base.count(x->x==n,p[i])) for n ∈ Base.unique(p[i])])

    for j ∈ Base.eachindex(Base.unique(p[i]))
        o = o ÷ Base.factorial(d[Base.unique(p[i])[j]])
    end

    return o

end


# ============================== triangle_coefficient(a, b, c) =============================

@doc raw"""
    triangle_coefficient(a::Real, b::Real, c::Real)
Triangle coefficient for a triangle of sides `a`, `b` and `c`.
#### Example:
```
julia> triangle_coefficient(3, 4, 5)
1//180180
julia> triangle_coefficient(1//2, 1, 1.5)
1//12
```
"""
function triangle_coefficient(a::Real, b::Real, c::Real)

    (a,b,c) = promote(a,b,c)

    isinteger(a + b + c) || return 0

    A = Int(a + b - c)
    B = Int(b + c - a)
    C = Int(c + a - b)

    A = A ≥ 0 ? bigfactorial(A) : return 0
    B = B ≥ 0 ? bigfactorial(B) : return 0
    C = C ≥ 0 ? bigfactorial(C) : return 0

    num = A * B * C
    den = bigfactorial(Int(a+b+c+1))

    return num//den

end

# ============================ istriangle(a, b, c) =============================

@doc raw"""
    istriangle(a::Real, b::Real, c::Real)
Triangle condition for a triangle of sides `a`, `b` and `c`.
#### Example:
```
julia> istriangle(3, 4, 5)
true
julia> istriangle(1//2, 1, 1.5)
true
```
"""
function istriangle(a::Real, b::Real, c::Real)

    Δ = triangle_coefficient(a,b,c)

    valid = Δ > 0 ? true : false

    return valid

end

# ...................... texp(x, p) .........................................

function _texp_int(x, p::Int)

    o = y = typeof(x)(1)

    x ≠ 0 || return o

    for n=1:p
        y *= x//n
        o += y
    end

    return o

end

function _texp_real(x, p::Int)

    o = y = typeof(x)(1)

    x ≠ 0.0 || return o

    for n=1:p
        y *= x/n
        o += y
    end

    return o

end

@doc raw"""
    texp(x::T, a::T, p::Int) where T <: Real
Taylor expansion of exp(x) about ``x = a`` up to order p.
```math
    \mathsf{texp}(x,a,p) = 1 + (x-a) + \frac{1}{2}(x-a)^2 + ⋯ + \frac{1}{p!}(x-a)^p.
```
### Examples:
```
julia> p = 5;
julia> texp(1.0, 0.0, 5)
2.7166666666666663
julia> texp(1, 0, 5)
163//60
```
"""
function texp(x::T, a::T, p::Int) where T <: Real

    x = x - a

    V = typeof(x)

    return  V <: Rational ? _texp_int(x, p) : V <: Integer ? _texp_int(x, p) : _texp_real(x, p)

end