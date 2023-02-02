# ==============================================================================
#                          Faulhaber.jl
# ==============================================================================

@doc raw"""
    faulhaber_polynom(p [, T=Int])

Vector representation of the Faulhaber polynomial of degree ``p``,
```math
    F(n,p)=\frac{1}{p}\sum_{j=1}^{p}{\binom {p}{p-j}}B_{p-j}n^{j}.
```
``F(n,p)=`` `polynomial(c,n)`, where ``c=[c_0,⋯\ c_p]`` is the coefficient vector, with
```math
    c_0=0,\ c_j=\frac{1}{p}{\binom {p}{p-j}}B_{p-j},
```
with ``j∈\{ 1,⋯\ p\}``. The ``B_0,⋯\ B_{p-1}`` are Bernoulli numbers
(but with ``B_1=+\frac{1}{2}`` rather than ``-\frac{1}{2}``).
### Example:
```
faulhaber_polynom(6)
7-element Vector{Rational{Int64}}:
  0//1
  0//1
 -1//12
  0//1
  5//12
  1//2
  1//6
```
"""
function faulhaber_polynom(k::Int; T=Int)

    k < 1 && return 0
    k > 1 || return 1 // 1

    P = CamMath.pascal_triangle(k)[end][1:end-1]
    B = CamMath.bernoulliB_array(k - 1)
    B[2] = -B[2]  # was bernoulliB_array(k-1; T)

    F = (B .* P) // k

    F = Base.append!(F, 0 // 1)   # add polynomial constant (zero in this case)

    return Base.reverse(F)     # reverse to standard order

end

# =================================== faulhaber_summation(n,p;T) ===============

@doc raw"""
    faulhaber_summation(n, p [, T=Int])

Sum of powers of natural numbers ``1,⋯\ n``,
```math
    FS(n,p)=\sum_{k=1}^{n}k^{p}=F(n,p+1).
```
where ``F(n,p)`` is the Faulhamer polynomial of degree ``p``.
### Examples:
```
faulhaber_summation(5,1)
 15

faulhaber_summation(3,60; T=BigInt)
  42391158276369125018901280178
```
"""
function faulhaber_summation(n::Int, p::Int)   # short argument: better performance

    n ≠ 0 || return 0

    F = CamMath.faulhaber_polynom(p + 1)
    o = 0
    for k = 1:p+1
        for i = 1:k
            F[k+1] *= n # avoid n^k in o = Base.sum([F[k+1]*n^k for k=1:p+1])
        end
        o += F[k+1]
    end

    Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")

    return Base.numerator(o)

end
function faulhaber_summation(n::Int, p::Int; T=Int)

    n ≠ 0 || return nothing

    F = CamMath.faulhaber_polynom(p + 1; T)
    o = 0
    for k = 1:p+1
        for i = 1:k
            F[k+1] *= n # avoid n^k in o = Base.sum([F[k+1]*n^k for k=1:p+1])
        end
        o += F[k+1]
    end

    Base.denominator(o) == 1 || error("Error: Faulhaber sum failed")

    return Base.numerator(o)

end