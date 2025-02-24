# SPDX-License-Identifier: MIT

# Copyright (c) 2023 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# ==============================================================================
#                            Bernoulli.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#               bernoulliB(n::Integer [; arr=false [, msg=true]])
# ------------------------------------------------------------------------------

function _bernoulli_BigInt(n, o)

    nul = big(0)
    one = big(1)
    l = length(o)
    o = collect(o)       # [o[n] for n = 1:l]    # transform tuple to array

    for k = l+1:n+1
        a = nul
        if Base.isodd(k)
            b = one
            for j = 1:k-1
                a -= o[j] * b
                b *= (k + 1 - j)
                b ÷= j                 # binomial coefficients are integers
            end
        end
        Base.push!(o, a // big(k))
    end

    return o

end

# ..............................................................................
@doc raw"""
    bernoulliB(n::Integer [; arr=false [, msg=true]])

Bernoulli numbers of index `n` are defined by the recurrence relation
```math
    B_n = - \frac{1}{n+1}\sum_{k=0}^{n-1}\frac{(n+1)!}{k!(n+1-k)}B_k,
```
with ``B_0=1`` and ``B_1=-1/2``. Including ``B_0`` results in the *even index 
convention* ``(B_{2n+1}=0`` for ``n>1)``.

- `arr` : output in array format

- `msg` : integer-overflow protection (IOP) - warning on activation 
#### Examples:
```
julia> o = [bernoulliB(n) for n=0:5]; println(o)
Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]

julia> o = bernoulliB(5; arr=true); println(o)
Rational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]

julia> o = bernoulliB(big(5); arr=true); println(o)
Rational{BigInt}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]

julia> bernoulliB(60)
IOP capture: bernoulliB(60) converted to Rational{BigInt}
-1215233140483755572040304994079820246041491//56786730

julia> n = 60;
julia> bernoulliB(n; msg=false) == bernoulliB(n; msg=false, arr=true)[end]             
true
```
"""
function bernoulliB(n::Integer; arr=false, msg=true)

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

    nc = 35
    no = 86

    T = Type_IOP(n, nc; nam="bernoulliB", msg)

    n = convert(Int, n)
    n ≥ 0 || throw(DomainError(n))

    if arr
        if n ≤ nc
            return (num[1:1+n] .// den[1:1+n])
        elseif n ≤ no
            return (num[1:1+n] .// den[1:1+n])
        else
            o = num .// den
            return _bernoulli_BigInt(n, o)
        end
    else
        if n ≤ nc
            return Rational{T}(num[1+n], den[1+n])
        elseif n ≤ no
            iseven(n) || return big(0) // big(1)
            return Rational{BigInt}(num[1+n], den[1+n])
        else
            iseven(n) || return big(0) // big(1)
            o = num .// den
            return _bernoulli_BigInt(n, o)[end]
        end
    end

end

