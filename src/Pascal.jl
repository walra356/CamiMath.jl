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
#                               Pascal.jl
#                           Jook Walraven - 14-2-2023
# ==============================================================================

# ------------------------------------------------------------------------------
#                           pascal_next(row)
# ------------------------------------------------------------------------------

@doc raw"""
    pascal_next(row)

Next `row` of binomial coefficients of the Pascal triangle. 
#### Example:
```
julia> pascal_next([1, 4, 6, 4, 1])
6-element Vector{Int64}:
  1
  5
 10
 10
  5
  1
```
"""
function pascal_next(a)

    n = Base.length(a) + 1
    o = Base.ones(eltype(a), n)
    n ≠ 2 || return o

    for k = 1:n÷2
        o[k+1] = a[k+1] + a[k]
        o[n-k] = o[k+1]
    end

    return o

end

# ------------------------------------------------------------------------------
#           pascal_triangle(n::Integer [[; arr=false], msg=true])
# ------------------------------------------------------------------------------

@doc raw"""
    pascal_triangle(n::Integer [; arr=false [, msg=true]])

Row `n` of Pascal triangle, ``r_n = [\binom{n}{1},\cdots\ \binom{n}{n}]``

- `arr` : output full Pascal triangle

- `msg` : integer-overflow protection (IOP) - warning on activation 
(for `n > 10000`)
#### Examples:
```
julia> [pascal_triangle(n) for n=0:5]
6-element Vector{Vector{Int64}}:
 [1]
 [1, 1]
 [1, 2, 1]
 [1, 3, 3, 1]
 [1, 4, 6, 4, 1]
 [1, 5, 10, 10, 5, 1]

julia> pascal_triangle(5; arr=true)
5-element Vector{Vector{Int64}}:
 [1, 1]
 [1, 2, 1]
 [1, 3, 3, 1]
 [1, 4, 6, 4, 1]
 [1, 5, 10, 10, 5, 1]
```
"""
function pascal_triangle(n::Integer; arr=false, msg=true)

    o = (
        (1, 1),
        (1, 2, 1),
        (1, 3, 3, 1),
        (1, 4, 6, 4, 1),
        (1, 5, 10, 10, 5, 1),
        (1, 6, 15, 20, 15, 6, 1),
        (1, 7, 21, 35, 35, 21, 7, 1),
        (1, 8, 28, 56, 70, 56, 28, 8, 1),
        (1, 9, 36, 84, 126, 126, 84, 36, 9, 1),
        (1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1),
        (1, 11, 55, 165, 330, 462, 462, 330, 165, 55, 11, 1),
        (1, 12, 66, 220, 495, 792, 924, 792, 495, 220, 66, 12, 1),
        (1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1),
        (1, 14, 91, 364, 1001, 2002, 3003, 3432, 3003, 2002, 1001, 364,
            91, 14, 1),
        (1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365,
            455, 105, 15, 1),
        (1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368,
            1820, 560, 120, 16, 1),
        (1, 17, 136, 680, 2380, 6188, 12376, 19448, 24310, 24310, 19448,
            12376, 6188, 2380, 680, 136, 17, 1),
        (1, 18, 153, 816, 3060, 8568, 18564, 31824, 43758, 48620, 43758,
            31824, 18564, 8568, 3060, 816, 153, 18, 1),
        (1, 19, 171, 969, 3876, 11628, 27132, 50388, 75582, 92378, 92378,
            75582, 50388, 27132, 11628, 3876, 969, 171, 19, 1),
        (1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756,
            167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1),
        (1, 21, 210, 1330, 5985, 20349, 54264, 116280, 203490, 293930,
            352716, 352716, 293930, 203490, 116280, 54264, 20349, 5985, 1330,
            210, 21, 1),
        (1, 22, 231, 1540, 7315, 26334, 74613, 170544, 319770, 497420, 646646,
            705432, 646646, 497420, 319770, 170544, 74613, 26334, 7315, 1540,
            231, 22, 1),
        (1, 23, 253, 1771, 8855, 33649, 100947, 245157, 490314, 817190,
            1144066, 1352078, 1352078, 1144066, 817190, 490314, 245157,
            100947, 33649, 8855, 1771, 253, 23, 1),
        (1, 24, 276, 2024, 10626, 42504, 134596, 346104, 735471, 1307504,
            1961256, 2496144, 2704156, 2496144, 1961256, 1307504, 735471,
            346104, 134596, 42504, 10626, 2024, 276, 24, 1),
        (1, 25, 300, 2300, 12650, 53130, 177100, 480700, 1081575, 2042975,
            3268760, 4457400, 5200300, 5200300, 4457400, 3268760, 2042975,
            1081575, 480700, 177100, 53130, 12650, 2300, 300, 25, 1))

    no = 25

    T = Type_IOP(n, 10000, "pascal_triangle($n)"; msg)
    
    n = convert(Int, n)
    n ≠ 0 || return arr ? [[T(1)]] : [T(1)]
    n > 0 || throw(DomainError(n))

    if arr
        if n ≤ no
            return collect(collect.(T, o))[1:n]
        else
            o = collect(collect.(T, o))
            for i = 1:n-no
                a = pascal_next(o[end])
                o = push!(o, a)
            end
            return o
        end
    else
        if n ≤ no
            return collect(T, o[n])
        else
            o = collect(T, o[no])
            for i = 1:n-no
                o = pascal_next(o)
            end
            return o
        end
    end

end