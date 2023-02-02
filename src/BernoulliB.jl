# SPDX-License-Identifier: MIT

# =============================================================================
#                            BernoulliB.jl
# =============================================================================

global glBn_Int = [1 // 1, -1 // 2, 1 // 6, 0 // 1, -1 // 30, 0 // 1, 1 // 42,
    0 // 1, -1 // 30, 0 // 1, 5 // 66, 0 // 1, -691 // 2730, 0 // 1,
    7 // 6, 0 // 1, -3617 // 510, 0 // 1, 43867 // 798, 0 // 1,
    -174611 // 330, 0 // 1, 854513 // 138, 0 // 1, -236364091 // 2730,
    0 // 1, 8553103 // 6, 0 // 1, -23749461029 // 870, 0 // 1,
    8615841276005 // 14322, 0 // 1, -7709321041217 // 510, 0 // 1,
    2577687858367 // 6, 0 // 1]

# ..............................................................................

global glBn_BigInt = convert(Vector{Rational{BigInt}}, glBn_Int)

# ..............................................................................

function _bn_BigInt(n::Int, nc::Int)

    nul = big(0)
    one = big(1)

    o = glBn_BigInt[1:1+nc]
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

Bernoulli numbers are defined by the recurrence relation
```math
    B_n = - \frac{1}{n+1}\sum_{k=0}^{n-1}\frac{(n+1)!}{k!(n+1-k)}B_k,
```
with ``B_0=1`` and ``B_1=-1/2``. Starting at ``B_0`` is called the *even index 
convention* since ``B_{2n+1}=0`` for ``n>1``.

Integer-overflow protection: for `n > 35` the output is autoconverted to 
`Rational{BigInt}`. By default the capture message is activated: 
"Warning: bernoulliB autoconverted to Rational{BigInt}". 
### Examples:
```
julia> bernoulliB(60)
"Warning: bernoulliB autoconverted to Rational{BigInt}"
-1215233140483755572040304994079820246041491//56786730

julia> bernoulliB_array(10; println(o)
Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30, 0//1]

julia> n = 60;
julia> bernoulliB(n) == bernoulliB_array(n)[end]             
true
```
"""
function bernoulliB(n::T; msg=true) where {T<:Integer}

    o = CamMath.bernoulliB_array(n; msg)[end]

    return o

end

# ..............................................................................

@doc raw"""
    bernoulliB_array(nmax::T [; msg=true]) where {T<:Integer}

Bernoulli number array, ``[B_0,⋯\ B_{nmax}]``. This is the *even index 
convention* (odd-indexed numbers vanish - except B[1]=-1/2). The array is 
calculated by repetative use of the recurrence relation.
```math
    B_n = - \frac{1}{n+1}\sum_{k=0}^{n-1}\frac{(n+1)!}{k!(n+1-k)}B_k.
```
Special numbers: ``B_0=1,\ B_1=-1/2,\ B_{2n+1}=0\ (\rm{for}\ n>1)``. Starting 
at ``B_0`` is called the *even index convention*.

Integer-overflow protection (IOP): for `n > 35` the output is autoconverted to 
`Rational{BigInt}`. By default the capture message is activated: 
"Integer-overflow protection: bernoulliB converted to Rational{BigInt}"
### Examples:
```
julia> "Integer-overflow protection: bernoulliB converted to Rational{BigInt}"
"Warning: bernoulliB autoconverted to Rational{BigInt}"
-1215233140483755572040304994079820246041491//56786730

julia> o = bernoulliB_array(8); println(o)
Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30]

julia> n = 60;
julia> bernoulliB(n) == bernoulliB_array(n)[end]             
true
```
"""
function bernoulliB_array(nmax::T; msg=true) where {T<:Integer}

    str = "Warning: IOP - bernoulliB converted to Rational{BigInt}"

    n = Int(nmax)
    nc = 35

    if n ≤ nc
        o = T == Int ? glBn_Int[1:1+n] : glBn_BigInt[1:1+n]
    else
        o = _bn_BigInt(n, nc)
        msg && T == Int && println(str)
    end

    return o

end
