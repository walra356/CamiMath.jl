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
function _bernoulli_BigInt(n::Int, o)

    nul = big(0)
    one = big(1)
    no = length(o)
    o = [o[n] for n = 1:no]    # transform tuple to array

    for m = no+1:n+1
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

Integer overflow protection (IOP): for `n > 35` the output is converted to 
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
function bernoulliB(n::Integer; msg=true)

    num = (
        1, -1, 1, 0, -1, 0, 1, 0, -1, 0, 5, 0, -691, 0, 7, 0, -3617, 0, 43867,
        0, -174611, 0, 854513, 0, -236364091, 0, 8553103, 0, -23749461029,
        0, 8615841276005, 0, -7709321041217, 0, 2577687858367,
        0, -26315271553053477373, 0, 2929993913841559,
        0, -261082718496449122051, 0, 1520097643918070802691,
        0, -27833269579301024235023, 0, 596451111593912163277961,
        0, -5609403368997817686249127547, 0, 495057205241079648212477525,
        0, -801165718135489957347924991853, 0, 29149963634884862421418123812691,
        0, -2479392929313226753685415739663229,
        0, 84483613348880041862046775994036021,
        0, -1215233140483755572040304994079820246041491,
        0, 12300585434086858541953039857403386151,
        0, -106783830147866529886385444979142647942017,
        0, 1472600022126335654051619428551932342241899101,
        0, -78773130858718728141909149208474606244347001,
        0, 1505381347333367003803076567377857208511438160235,
        0, -5827954961669944110438277244641067365282488301844260429,
        0, 34152417289221168014330073731472635186688307783087,
        0, -24655088825935372707687196040585199904365267828865801,
        0, 414846365575400828295179035549542073492199375372400483487,
        0, -4603784299479457646935574969019046849794257872751288919656867,
        0, 1677014149185145836823154509786269900207736027570253414881613,
        0, -2024576195935290360231131160111731009989917391198090877281083932477,
        0, 660714619417678653573847847426261496277830686653388931761996983
    )

    den = (
        1, 2, 6, 1, 30, 1, 42, 1, 30, 1, 66, 1, 2730, 1, 6, 1, 510, 1, 798,
        1, 330, 1, 138, 1, 2730, 1, 6, 1, 870, 1, 14322, 1, 510, 1, 6,
        1, 1919190, 1, 6, 1, 13530, 1, 1806, 1, 690, 1, 282, 1, 46410, 1, 66,
        1, 1590, 1, 798, 1, 870, 1, 354, 1, 56786730, 1, 6, 1, 510, 1, 64722,
        1, 30, 1, 4686, 1, 140100870, 1, 6, 1, 30, 1, 3318, 1, 230010, 1, 498,
        1, 3404310, 1, 6
    )

    str = "IOP capture: bernoulliB(n) converted to Rational{BigInt}"

    T = n > 35 ? BigInt : typeof(n)

    n′ = 1 + n

    iseven(n) || Rational{typeof(n)}(0, 1)

    if n < 0
        throw(DomainError(n))
    elseif n ≤ 35
        return Rational{T}(num[n′], den[n′])
    elseif n ≤ 86
        msg && typeof(n) == Int && println(str)
        return Rational{T}(num[n′], den[n′])
    else
        msg && typeof(n) == Int && println(str)
        o = Rational{T}.(num, den)
        return _bernoulli_BigInt(n, o)[end]
    end

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
        o = _bernoulli_BigInt(1+n, nc)
        msg && T == Int && println(str)
    end

    return o

end


function bernoulliB_array1(nmax::Integer; msg=true)

     o = ( 
        1 // 1, -1 // 2, 1 // 6, 0 // 1, -1 // 30, 0 // 1,
        1 // 42, 0 // 1, -1 // 30, 0 // 1, 5 // 66, 0 // 1, -691 // 2730,
        0 // 1, 7 // 6, 0 // 1, -3617 // 510, 0 // 1, 43867 // 798, 0 // 1,
        -174611 // 330, 0 // 1, 854513 // 138, 0 // 1, -236364091 // 2730,
        0 // 1, 8553103 // 6, 0 // 1, -23749461029 // 870, 0 // 1,
        8615841276005 // 14322, 0 // 1, -7709321041217 // 510, 0 // 1,
        2577687858367 // 6, 0 // 1
     )

    str = "IOP capture: bernoulliB_array(nmax) converted to Rational{BigInt}"

    T = nmax > 35 ? BigInt : typeof(nmax)

    n = convert(Int, nmax)

    if n < 0
        throw(DomainError(n))
    elseif n ≤ 35
        return @view o[1:1+n]
    else
        msg && T == Int && println(str)
        return _bernoulli_BigInt(1 + n, o)
    end

end

function _bernoulli_BigInt1(n::Int, nc::Int)

    nul = big(0)
    one = big(1)

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