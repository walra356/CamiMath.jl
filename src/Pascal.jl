# ==================================== pascal_triangle(nmax)  ============

@doc raw"""
    pascal_triangle(nmax [, T=Int])

Pascal triangle of binomial coefficients ``\binom{n}{k}`` for ``n=0,\ 1,⋯\ nmax``
### Example:
```
pascal_triangle(5)
6-element Vector{Vector{Int64}}:
 [1]
 [1, 1]
 [1, 2, 1]
 [1, 3, 3, 1]
 [1, 4, 6, 4, 1]
 [1, 5, 10, 10, 5, 1]
```
"""
function pascal_triangle(nmax::T) where {T<:Integer}

    nmax < 0 && error("Error: nmax must be a non-negative integer")
    nmax > T(10000) && error("Error: integer overflow")

    o = [Base.ones(T, n + 1) for n = 0:nmax]

    for n = 2:nmax
        for k = 1:n÷2
            o[n+1][k+1] = o[n][k+1] + o[n][k]
            o[n+1][n+1-k] = o[n+1][k+1]
        end
    end

    return o

end

# ==================================== pascal_next(nmax)  ======================

@doc raw"""
    pascal_next(nmax)

Next row of Pascal triangle
### Example:
```
a = [1, 4, 6, 4, 1]
pascal_next(a)
 [1, 5, 10, 10, 5, 1]
```
"""
function pascal_next(a::Vector{Int})

    n = Base.length(a) + 1
    o = Base.ones(Int,n)

    for k=1:n÷2
        o[k+1] = a[k+1] + a[k]
        o[n-k] = o[k+1]
    end

    return o

end