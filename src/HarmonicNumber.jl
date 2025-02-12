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
#                          harmonicNumber.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#        harmonicNumber(n::Integer [, p::Int [; arr=false [, msg=true]]])
# ------------------------------------------------------------------------------


function _harmonicNumbers(n::Int, p::Int)

    one = big(1)

    o = Rational{BigInt}[one//one]

    for k = 2:n
        push!(o, o[end] + one // big(k)^p)
    end

    return o

end
# ..............................................................................
function _harmonicNumbers_next(n::Int, nc::Int, p::Int, o)

    one = big(1)

    for k = nc+1:n
        push!(o, o[end] + one // big(k)^p)
    end

    return o

end

@doc raw"""
    harmonicNumber(n::Integer [, p=1 [; arr=false [, msg=true]]])

Sum of the ``p^{th}`` power of reciprocals of the first ``n`` positive integers,
```math
    H_{n,p}=\sum_{k=1}^{n}\frac{1}{k^p}.
```
- `arr` : output in array format

- `msg` : integer-overflow protection (IOP) - warning on activation 
#### Examples:
```
julia> o = [harmonicNumber(n) for n=1:8]; println(o)
Rational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280]

julia> harmonicNumber(8; arr=true)
(1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280)

julia> harmonicNumber(42)
12309312989335019//2844937529085600

julia> harmonicNumber(43)
IOP capture: harmonicNumber(43, 1) converted to Rational{BigInt}
532145396070491417//122332313750680800

julia> harmonicNumber(12) == harmonicNumber(12, 1)
true

julia> harmonicNumber(12, -3) == faulhaber_summation(12, 3)
true

julia> o = [harmonicNumber(i, 5) for i=1:4]; println(o)
Rational{Int64}[1//1, 33//32, 8051//7776, 257875//248832]

julia> o = harmonicNumber(4, 5; arr=true); println(o)
(1//1, 33//32, 8051//7776, 257875//248832)
```
"""
function harmonicNumber(n::Integer, p=1; arr=false, msg=true)

    n1 = (219060189739591200, 328590284609386800, 401610347855917200,
        456375395290815000, 500187433238733240, 536697464861998440,
        567991777681940040, 595374301399388940, 619714322481565740,
        641620341455524860, 661534904159124060, 679789919970756660,
        696640703796879060, 712287860206849860, 726891872856155940,
        740583134714880390, 753469028228973990, 765639038770062390,
        777168522440567190, 788121531927546750, 798552969534193950,
        808510250885993550, 818034606961627950, 827162114867444250,
        835924522457027898, 844349914370089098, 852463254730814698,
        860286832935800098, 867840632581992898, 875142638906645938,
        882209096640181138, 889054727569543363, 895692915137409763,
        902135861894456563, 908394724458444883, 914479729728989083,
        920400275397626683, 926165017232879083, 931781945174919883,
        937258449918409663, 942601381375472863, 947817100178796463)

    n2 = (54192375991353600, 67740469989192000, 73761845099342400,
        77148868598802000, 79316563638456144, 80821907415993744,
        81927874272960144, 82774630147825044, 83443671826730644,
        83985595586644180, 84433466462605780, 84809802406990180,
        85130467353684580, 85406959067926180, 85647814072332196,
        85859503041048421, 86047019912990821, 86214280332717221,
        86364397717734821, 86499878657713205, 86622763864042805,
        86734731583033205)

    n3 = (374368864117248000, 421164972131904000, 435030485617728000,
        440879999119560000, 443874950032497984, 445608139218225984,
        446699593632561984, 447430782820290984, 447944320356802984,
        448318689220920232, 448599958089528232, 448816606737744232,
        448987006766928232, 449123438568720232, 449234362676606824,
        449325761325072949)

    n4 = (590436101122560000, 627338357442720000, 634627692024480000,
        636934083044490000, 637878780806286096, 638334364217646096,
        638580276796206096, 638724426234956721, 638814418019916721,
        638873461630028977, 638913789210188977, 638942263173398977)

    n5 = (101625502003200000, 104801298940800000, 105219510883200000,
        105318754537500000, 105351274698141024, 105364343821341024,
        105370390438941024, 105373491803137899, 105375212839937899,
        105376229094957931)
    n6 = (351298031616000000, 356787063360000000, 357268953664000000,
        357354719785000000, 357377202859023424, 357384732395023424,
        357387718379023424, 357389058474664049)

    n7 = (2305393332480000000, 2323404217890000000, 2324458352930000000,
        2324599062972265625, 2324628572006921369, 2324636807436921369,
        2324639606796921369)

    n8 = (167961600000000, 168617700000000, 168643300000000, 168645862890625,
        168646292872321, 168646392872321)

    n9 = (10077696000000000, 10097379000000000, 10097891000000000,
        10097929443359375, 10097934603139727, 10097935603139727)

    n0 = (604661760000000000, 605252250000000000, 605262490000000000,
        605263066650390625, 605263128567754849, 605263138567754849)


    D = (219060189739591200, 54192375991353600, 374368864117248000,
        590436101122560000, 101625502003200000, 351298031616000000,
        2305393332480000000, 167961600000000, 10077696000000000,
        604661760000000000)

    no = (42, 22, 16, 12, 10, 8, 7, 6, 6, 6)     # no(p) == 6 for p=10

    n ≠ 0 || return n
    p ≠ 0 || return n

    if p > 0 # -----------------------------------------------------------------
        nc = p < 11 ? no[p] : p < 18 ? 4 : p < 25 ? 3 : 0

        n ≥ 0 || throw(DomainError(n))
        
        T = Type_IOP(n, nc, p; nam="harmonicNumber", msg)

        n = convert(Int, n)

        if arr # ...............................................................
            if n ≤ nc
                o = p == 1 ? (n1.//T(D[p]))[1:n] :
                    p == 2 ? (n2.//T(D[p]))[1:n] :
                    p == 3 ? (n3.//T(D[p]))[1:n] :
                    p == 4 ? (n4.//T(D[p]))[1:n] :
                    p == 5 ? (n5.//T(D[p]))[1:n] :
                    p == 6 ? (n6.//T(D[p]))[1:n] :
                    p == 7 ? (n7.//T(D[p]))[1:n] :
                    p == 8 ? (n8.//T(D[p]))[1:n] :
                    p == 9 ? (n9.//T(D[p]))[1:n] :
                    p == 10 ? (n0.//T(D[p]))[1:n] :
                    _harmonicNumbers(n, p)
                return o
            else
                o = p == 1 ? Rational{T}[n1[i] // D[1] for i = 1:nc] :
                    p == 2 ? Rational{T}[n2[i] // D[2] for i = 1:nc] :
                    p == 3 ? Rational{T}[n3[i] // D[3] for i = 1:nc] :
                    p == 4 ? Rational{T}[n4[i] // D[4] for i = 1:nc] :
                    p == 5 ? Rational{T}[n5[i] // D[5] for i = 1:nc] :
                    p == 6 ? Rational{T}[n6[i] // D[6] for i = 1:nc] :
                    p == 7 ? Rational{T}[n7[i] // D[7] for i = 1:nc] :
                    p == 8 ? Rational{T}[n8[i] // D[8] for i = 1:nc] :
                    p == 9 ? Rational{T}[n9[i] // D[9] for i = 1:nc] :
                    p == 10 ? Rational{T}[n0[i] // D[10] for i = 1:nc] :
                    _harmonicNumbers(nc, p)
                return _harmonicNumbers_next(n, nc, p, o)
            end
        else #..................................................................
            if n ≤ nc
                o = p == 1 ? (n1[n] // T(D[p])) :
                    p == 2 ? (n2[n] // T(D[p])) :
                    p == 3 ? (n3[n] // T(D[p])) :
                    p == 4 ? (n4[n] // T(D[p])) :
                    p == 5 ? (n5[n] // T(D[p])) :
                    p == 6 ? (n6[n] // T(D[p])) :
                    p == 7 ? (n7[n] // T(D[p])) :
                    p == 8 ? (n8[n] // T(D[p])) :
                    p == 9 ? (n9[n] // T(D[p])) :
                    p == 10 ? (n0[n] // T(D[p])) :
                    _harmonicNumbers(n, p)[end]
                return o
            else
                o = p == 1 ? Rational{T}(n1[end], D[1]) :
                    p == 2 ? Rational{T}(n2[end], D[2]) :
                    p == 3 ? Rational{T}(n3[end], D[3]) :
                    p == 4 ? Rational{T}(n4[end], D[4]) :
                    p == 5 ? Rational{T}(n5[end], D[5]) :
                    p == 6 ? Rational{T}(n6[end], D[6]) :
                    p == 7 ? Rational{T}(n7[end], D[7]) :
                    p == 8 ? Rational{T}(n8[end], D[8]) :
                    p == 9 ? Rational{T}(n9[end], D[9]) :
                    p == 10 ? Rational{T}(n0[end], D[10]) :
                    _harmonicNumbers(nc, p)[end]
                return _harmonicNumbers_next(n, nc, p, [o])[end]
            end
        end # ..................................................................
    else # ---------------------------------------------------------------------
        return CamiMath.faulhaber_summation(n, -p; msg)
    end  # ---------------------------------------------------------------------

end