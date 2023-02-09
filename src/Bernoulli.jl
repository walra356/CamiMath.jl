# SPDX-License-Identifier: MIT

# ==============================================================================
#                            Bernoulli.jl
# ==============================================================================

global gl_bernou_Int = [
    1 // 1, -1 // 2, 1 // 6, 0 // 1, -1 // 30, 0 // 1,
    1 // 42, 0 // 1, -1 // 30, 0 // 1, 5 // 66, 0 // 1, -691 // 2730,
    0 // 1, 7 // 6, 0 // 1, -3617 // 510, 0 // 1, 43867 // 798, 0 // 1,
    -174611 // 330, 0 // 1, 854513 // 138, 0 // 1, -236364091 // 2730,
    0 // 1, 8553103 // 6, 0 // 1, -23749461029 // 870, 0 // 1,
    8615841276005 // 14322, 0 // 1, -7709321041217 // 510, 0 // 1,
    2577687858367 // 6, 0 // 1
]

# ..............................................................................

global gl_bernou_BigInt = convert(Vector{Rational{BigInt}}, gl_bernou_Int)

# ..............................................................................

function _bernoulli_BigInt(n::Int, nc::Int)

    nul = big(0)
    one = big(1)

    o = gl_bernou_BigInt[1:1+nc]
    for m = nc+2:n+1
        a = nul
        if Base.isodd(m)
            b = one
            for j = 1:m-1
                a -= o[j] * b
                b *= (m + 1 - j)
                b ÷= j                     # binomial coefficients are integers
            end
        end
        Base.push!(o, a // big(m))
    end

    return o

end

# ..............................................................................

@doc raw"""
    bernoulliB(n::T [; msg=true]) where {T<:Integer}

Bernoulli numbers of index `n` are defined by the recurrence relation
```math
    B_n = - \frac{1}{n+1}\sum_{k=0}^{n-1}\frac{(n+1)!}{k!(n+1-k)}B_k,
```
with ``B_0=1`` and ``B_1=-1/2``. Including ``B_0`` results in the *even index 
convention* ``(B_{2n+1}=0`` for ``n>1)``.

Integer overflow protection (IOP): f0r `n > 35` the output is converted to 
Rational{BigInt}. By default the IOP capture message is activated.
### Examples:
```
julia> o = [bernoulliB(n) for n=0:5]; println(o)
Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]

julia> bernoulliB(60)
IOP capture: bernoulliB(60) converted to Rational{BigInt}
-1215233140483755572040304994079820246041491//56786730

julia> n = 60;
julia> bernoulliB(n) == bernoulliB_array(n)[end]             
true
```
"""
function bernoulliB(n::T; msg=true) where {T<:Integer}

    n ≠ 1 || return -T(1) // T(2)

    iseven(n) || return T(0) // T(1)

    o = CamiMath.bernoulliB_array(n; msg)[end]

    return o

end

# ..............................................................................

@doc raw"""
    bernoulliB_array(nmax::T [; msg=true]) where {T<:Integer}

Bernoulli number array ``[B_0,\cdots\ B_{nmax}]``, where `nmax` is the index of 
the highest Bernoulli number of the array (NB.: *not* the array length).
### Examples:
```
julia> o = bernoulliB_array(8); println(o)
Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30]

julia> o = bernoulliB_array(big(8)); println(o)
Rational{BigInt}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30]

julia> n = 60; msg = false;
julia>  bernoulliB(n; msg) == bernoulliB_array(n; msg)[end]            
true
```
"""
function bernoulliB_array(nmax::T; msg=true) where {T<:Integer}

    str = "IOP capture: bernoulliB($(nmax)) converted to Rational{BigInt}"

    n = Int(nmax)
    nc = 35

    if n ≤ nc
        o = T == Int ? gl_bernou_Int[1:1+n] : gl_bernou_BigInt[1:1+n]
    else
        o = _bernoulli_BigInt(n, nc)
        msg && T == Int && println(str)
    end

    return o

end
