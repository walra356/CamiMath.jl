# SPDX-License-Identifier: MIT

# author: Jook Walraven - 17-2-2023

# ==============================================================================
#                             Partitions.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#     canonical_partitions(n::Int, [[[m=0]; header=true,] reverse=true])
# ------------------------------------------------------------------------------

function _canonical_partition(n::Int, m::Int)

    o = Base.fill(m, Base.cld(n, m))           # init partition
    o[Base.cld(n, m)] = ((n % m) ≠ 0 ? n % m : m)
    # adjust last element of partition

    return o

end

@doc raw"""
    canonical_partitions(n::Int, [m=0 [; header=true [, reverse=true]]])

Canonical partition of `n` in parts of maximum size `m` (`m` = 0 for any size)

`header` : unit partition included in output
#### Examples:
```
julia> canonical_partitions(6; header=true, reverse=false)
6-element Vector{Vector{Int64}}:
 [6]
 [5, 1]
 [4, 2]
 [3, 3]
 [2, 2, 2]
 [1, 1, 1, 1, 1, 1]

julia> canonical_partitions(6; header=true)
6-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]
 [6]

julia> canonical_partitions(6)
6-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]
 [6]

julia> o = canonical_partitions(9, 2); println(o)
[2, 2, 2, 2, 1]

julia> o = canonical_partitions(9, 3); println(o)
[3, 3, 3]
```
"""
function canonical_partitions(n::Int, m=0, sense=rev; header=true)

    h = header ? n : n - 1

    if m < 0
        throw(DomainError(m))
    elseif m == 0
        if isreverse(sense)
            o = [_canonical_partition(n, m) for m = 1:h]
        else
            o = [_canonical_partition(n, m) for m = h:-1:1]
        end
    elseif 0 < m ≤ n
        o = _canonical_partition(n, m)
    else
        throw(DomainError(m))
    end

    return o

end
function canonical_partitions1(n::Int, m=0; header=true, reverse=true)

    h = header ? n : n - 1

    if m < 0
        throw(DomainError(m))
    elseif m == 0
        if reverse
            o = [_canonical_partition(n, m) for m = 1:h]
        else
            o = [_canonical_partition(n, m) for m = h:-1:1]
        end
    elseif 0 < m ≤ n
        o = _canonical_partition(n, m)
    else
        throw(DomainError(m))
    end

    return o

end

# ------------------------------------------------------------------------------
#      integer_partitions(n [[[,m]; transpose=false], count=false])
# ------------------------------------------------------------------------------

function _partition_count(n::Int, k::Int)

    o = (n < 0) | (k < 0) | (k > n) ? 0 :
        (k == n) | (k == 1) ? 1 :
        _partition_count(n - k, k) + _partition_count(n - 1, k - 1)

    return (o)

end

# ..............................................................................
function _partition(a::Vector{Int}, n::Int, i::Int, cp::Vector{Vector{Vector{Int}}})

    o = a[1:i-1]
    m = a[i] - 1                 # m: partition value
    ni = n - Base.sum(o)         # ni: sub-partition index at partition index i

    Base.append!(o, cp[ni][m])   # complete partition by appending it to a

    return o

end

# ..............................................................................
function _restricted_partitions(o::Vector{Int}, n::Int, np::Int, cp::Vector{Vector{Vector{Int}}})

    oo = [o]

    for p = 1:np-1
        i = Base.findlast(x -> x > 1, oo[p])
        Base.append!(oo, [_partition(oo[p], n, i, cp)])
    end

    return oo

end

# ..............................................................................
@doc raw"""
    integer_partitions(n [,m [; transpose=false [, count=false]]])

`default`                      : The integer partitions of `n`

`count`                        : The number of integer partitions of `n`

`transpose` = `false` (`m` > 0): partitions restricted to maximum part `m`
            = `true`  (`m` > 0): partitions restricted to maximum length `m``

definitions:

The integer partition of the positive integer `n` is a nonincreasing
sequence of positive integers p1, p2,... pk whose sum is `n`.
The elements of the sequence are called the parts of the partition.
#### Examples:
```
julia> integer_partitions(7)
15-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 1, 1, 1]
 [2, 2, 2, 1]
 [2, 2, 1, 1, 1]
 [2, 1, 1, 1, 1, 1]
 [3, 3, 1]
 [3, 2, 2]
 [3, 2, 1, 1]
 [3, 1, 1, 1, 1]
 [4, 3]
 [4, 2, 1]
 [4, 1, 1, 1]
 [5, 2]
 [5, 1, 1]
 [6, 1]
 [7]

julia> integer_partitions(7; count=true)
15

julia> integer_partitions(7, 4; count=true)
3

julia> integer_partitions(7, 4)
3-element Vector{Vector{Int64}}:
 [4, 3]
 [4, 2, 1]
 [4, 1, 1, 1]

julia> integer_partitions(7, 4; transpose=true)
3-element Vector{Vector{Int64}}:
 [2, 2, 2, 1]
 [3, 2, 1, 1]
 [4, 1, 1, 1]
```
"""
function integer_partitions(n::Int, m=0; transpose=false, count=false)

    cp = [canonical_partitions(m) for m = 1:n]
    pc = [_partition_count(n, m) for m = 1:n]
    oo = [ones(Int, n)]

    np = m > 0 ? pc[m] : sum(pc)

    if !count

        if m == 0
            o = [_restricted_partitions(cp[n][p], n, pc[p], cp) for p = 2:n]
            for p = 1:n-1
                append!(oo, o[p])
            end
        else
            oo = _restricted_partitions(cp[n][m], n, pc[m], cp)
        end

        if transpose
            for p = 1:np
                l = length(oo[p])
                s = max(oo[p][1], l)
                mat = zeros(Int, s, s)
                for j = 1:l
                    for i = 1:oo[p][j]
                        mat[i, j] = 1
                    end
                end
                oo[p] = [sum(mat[i, :]) for i = 1:oo[p][1]]
            end

        end

    end

    return count ? np : oo

end



