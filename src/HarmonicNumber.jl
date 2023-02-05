# SPDX-License-Identifier: MIT

# ==============================================================================
#                          harmonicNumber.jl
# ==============================================================================

global gl_harmon_Int = Vector{Rational{Int64}}[
    [1 // 1, 3 // 2, 11 // 6, 25 // 12, 137 // 60, 49 // 20, 363 // 140,
        761 // 280, 7129 // 2520, 7381 // 2520, 83711 // 27720, 86021 // 27720,
        1145993 // 360360, 1171733 // 360360, 1195757 // 360360,
        2436559 // 720720, 42142223 // 12252240, 14274301 // 4084080,
        275295799 // 77597520, 55835135 // 15519504, 18858053 // 5173168,
        19093197 // 5173168, 444316699 // 118982864, 1347822955 // 356948592,
        34052522467 // 8923714800, 34395742267 // 8923714800,
        312536252003 // 80313433200, 315404588903 // 80313433200,
        9227046511387 // 2329089562800, 9304682830147 // 2329089562800,
        290774257297357 // 72201776446800, 586061125622639 // 144403552893600,
        53676090078349 // 13127595717600, 54062195834749 // 13127595717600,
        54437269998109 // 13127595717600, 54801925434709 // 13127595717600,
        2040798836801833 // 485721041551200,
        2053580969474233 // 485721041551200,
        2066035355155033 // 485721041551200,
        2078178381193813 // 485721041551200,
        85691034670497533 // 19914562703599200,
        12309312989335019 // 2844937529085600,
        532145396070491417 // 122332313750680800,
        5884182435213075787 // 1345655451257488800,
        5914085889685464427 // 1345655451257488800,
        5943339269060627227 // 1345655451257488800],
    [1 // 1, 5 // 4, 49 // 36, 205 // 144, 5269 // 3600, 5369 // 3600,
        266681 // 176400, 1077749 // 705600, 9778141 // 6350400,
        1968329 // 1270080, 239437889 // 153679680, 240505109 // 153679680,
        40799043101 // 25971865920, 40931552621 // 25971865920,
        205234915681 // 129859329600, 822968714749 // 519437318400,
        238357395880861 // 150117385017600, 238820721143261 // 150117385017600,
        86364397717734821 // 54192375991353600,
        17299975731542641 // 10838475198270720,
        353562301485889 // 221193371393280, 354019312583809 // 221193371393280,
        187497409728228241 // 117011293467045120,
        187700554334941861 // 117011293467045120],
    [1 // 1, 9 // 8, 251 // 216, 2035 // 1728, 256103 // 216000,
        28567 // 24000, 9822481 // 8232000, 78708473 // 65856000,
        19148110939 // 16003008000, 19164113947 // 16003008000,
        25523438671457 // 21300003648000, 25535765062457 // 21300003648000,
        56123375845866029 // 46796108014656000,
        56140429821090029 // 46796108014656000,
        56154295334575853 // 46796108014656000,
        449325761325072949 // 374368864117248000],
    [1 // 1, 17 // 16, 1393 // 1296, 22369 // 20736, 14001361 // 12960000,
        14011361 // 12960000, 33654237761 // 31116960000,
        538589354801 // 497871360000, 43631884298881 // 40327580160000,
        43635917056897 // 40327580160000,
        638913789210188977 // 590436101122560000,
        638942263173398977 // 590436101122560000],
    [1 // 1, 33 // 32, 8051 // 7776, 257875 // 248832, 806108207 // 777600000,
        268736069 // 259200000,
        4516906311683 // 4356374400000,
        144545256245731 // 139403980800000,
        105375212839937899 // 101625502003200000,
        105376229094957931 // 101625502003200000],
    [1 // 1, 65 // 64, 47449 // 46656, 3037465 // 2985984,
        47463376609 // 46656000000, 47464376609 // 46656000000,
        5584183099672241 // 5489031744000000,
        357389058474664049 // 351298031616000000],
    [1 // 1, 129 // 128, 282251 // 279936, 36130315 // 35831808,
        2822716691183 // 2799360000000, 940908897061 // 933120000000,
        774879868932307123 // 768464444160000000],
    [1 // 1, 257 // 256, 1686433 // 1679616, 431733409 // 429981696,
        168646292872321 // 167961600000000,
        168646392872321 // 167961600000000],
    [1 // 1, 513 // 512, 10097891 // 10077696, 5170139875 // 5159780352,
        10097934603139727 // 10077696000000000,
        373997614931101 // 373248000000000],
    [1 // 1, 1025 // 1024, 60526249 // 60466176, 61978938025 // 61917364224,
        605263128567754849 // 604661760000000000,
        605263138567754849 // 604661760000000000]
]

# ..............................................................................
global gl_harmon_BigInt = convert.(Vector{Rational{BigInt}}, gl_harmon_Int)

# ..............................................................................
function _harmonicNumber_BigInt(n::Int, nc::Int)

    one = big(1)

    o = gl_harmon_BigInt[1][1:nc]
    for m = nc+1:n
        a = o[m-1] + one // big(m)
        Base.push!(o, a)
    end

    return o

end

# ..............................................................................
@doc raw"""
    harmonicNumber(n::T [; msg=true]) where {T<:Integer} 
   
Sum of the reciprocals of the first ``n`` natural numbers
```math
    H_n=\sum_{k=1}^{n}\frac{1}{k}.
```
On integer overflow the output is converted to Rational{BigInt}. 
By default the IOP capture message is activated: 
"Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}". 
### Examples:
```
julia> o = harmonicNumber_array(8); println(o)
Rational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280]

julia> o = [harmonicNumber(46; msg=true)]; println(o)
Rational{Int64}[5943339269060627227//1345655451257488800]

julia> o = [harmonicNumber(47; msg=true)]; println(o)
Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}
Rational{BigInt}[282057509927739620069//63245806209101973600]

julia> harmonicNumber(12) == harmonicNumber(12, 1)
true
```
"""
function harmonicNumber(n::T; msg=true) where {T<:Integer}

    o = CamMath.harmonicNumber_array(n; msg)[end]

    return o

end

# ==============================================================================
@doc raw"""
    harmonicNumber_array(nmax::T [; msg=true]) where {T<:Integer} 

Sum of the reciprocals of the first ``n`` natural numbers
```math
    H_n=\sum_{k=1}^{n}\frac{1}{k}.
```
On integer overflow the output is converted to 
Rational{BigInt}. By default the IOP capture message is activated: 
"Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}". 
### Examples:
```
julia> o = harmonicNumber_array(8); println(o)
Rational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280] 

julia> o = [harmonicNumber(46; msg=true)]; println(o)
Rational{Int64}[5943339269060627227//1345655451257488800]

julia> o = [harmonicNumber(47; msg=true)]; println(o)
Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}
Rational{BigInt}[282057509927739620069//63245806209101973600]

julia> harmonicNumber(12) == harmonicNumber(12, 1)
true
```
"""
function harmonicNumber_array(nmax::T; msg=true) where {T<:Integer}

    str1 = "Error: harmonic numbers defined for positive integers"
    str2 = "Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}"

    nmax > 0 || error(str1)

    n = Int(nmax)
    nc = 46

    if n ≤ nc
        o = T == Int ? gl_harmon_Int[1][1:n] : gl_harmon_BigInt[1][1:n]
    else
        o = _harmonicNumber_BigInt(n, nc)
        msg && T == Int && println(str2)
    end

    return o

end

# ==============================================================================
function _harmonicNumber_p_BigInt_header(n::Int, p::Int)

    nul = big(0)
    one = big(1)

    o = Rational{BigInt}[]
    if p > 10
        b = nul // one
        for m = 1:n
            a = one
            for i = 1:p
                a *= big(m)
            end
            b += one // a
            Base.push!(o, b)
        end
    else
        o = gl_harmon_BigInt[p]
    end

    return o

end
# ..............................................................................
function _harmonicNumber_p_BigInt(n::Int, nc::Int, p::Int)

    nul = big(0)
    one = big(1)

    c = CamMath._harmonicNumber_p_BigInt_header(nc, p)[1:nc]
    o = copy(c)

    #b = nul // one
    b = o[nc]
    for m = nc+1:n
        a = one
        for i = 1:p
            a *= big(m)
        end
        b += one // a
        Base.push!(o, b)
    end

    println(copy(o))
    return o

end
@doc raw"""
    harmonicNumber(n::T, p::Int [; msg=true]) where {T<:Integer}

Sum of the ``p^{th}`` power of reciprocals of the first ``n`` positive integers,
```math
    H_{n,p}=\sum_{k=1}^{n}\frac{1}{k^p}.
```
On integer overflow the output is converted to 
Rational{BigInt}. By default the IOP capture message is activated: 
"Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}". 
### Examples:
```
julia> o = [harmonicNumber(46,1; msg=true)]; println(o)
Rational{Int64}[5943339269060627227//1345655451257488800]

julia> o = [harmonicNumber(47,1; msg=true)]; println(o)
Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}
Rational{BigInt}[280682601097106968469//63245806209101973600]

julia> o = [harmonicNumber(47,1)]; println(o)
Rational{BigInt}[280682601097106968469//63245806209101973600]

harmonicNumber(12, -3) == faulhaber_summation(12, 3)
  true
```
"""
function harmonicNumber(n::T, p::Int; msg=true) where {T<:Integer}

    str = "Warning (IOP): harmonicNumber autoconverted to Rational{BigInt}"

    n ≠ 0 || return T(0)
    p ≠ 0 || return n
    p ≠ 1 || return harmonicNumber(n::T; msg) 

    if p > 0
        n = convert(Int, n)
        nc = p < 11 ? length(gl_harmon_Int[p]) : p < 18 ? 4 : p < 25 ? 3 : 0
        if (n ≤ nc) & (p < 11)
            println("A: n = $n nc = $(nc) p = $p")
            o = T == Int ? gl_harmon_Int[p][n] : gl_harmon_BigInt[p][n]
        else
            println("B: n = $n nc = $(nc) p = $p")
            o = _harmonicNumber_p_BigInt(n, nc, p)[n]
            msg && T == Int && println(str)
        end
    else
        o = faulhaber_summation(n, -p; msg)
    end

    return o

end