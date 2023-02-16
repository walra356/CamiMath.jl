# SPDX-License-Identifier: MIT

# author: Jook Walraven - 16-2-2023

# ==============================================================================
#                               Partitions.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#            fibonacci(n::Integer [[; arr=false], msg=true])
# ------------------------------------------------------------------------------

function _canonical_partition(n::Int, m::Int)

    o = Base.fill(m, Base.cld(n, m))           # init partition
    o[Base.cld(n, m)] = ((n % m) ≠ 0 ? n % m : m)
    # adjust last element of partition

    return o

end

@doc raw"""
    canonical_partitions(n [[; header=false], reverse=true])

"""
function canonical_partitions(n::Int, m=0; header=true, reverse=true)

    h = header ? n : n - 1

    if m == 0
        if reverse
            o = [_canonical_partition(n, m) for m = 1:h]
        else
            o = [_canonical_partition(n, m) for m = h:-1:1]
        end
    elseif 0 < m <= n
        o = _canonical_partition(n, m)
    else
        o = nothing
    end

    return o

end

function _partition_count(n::Int, k::Int)

    (n < 0) | (k < 0) | (k > n) ? 0 : (k == n) | (k == 1) ? 1 : _partition_count(n - k, k) + _partition_count(n - 1, k - 1)

end

function _partition(a::Vector{Int}, n::Int, i::Int, cp::Vector{Vector{Vector{Int}}})

    o = a[1:i-1]
    m = a[i] - 1                 # m: partition value
    ni = n - Base.sum(o)         # ni: sub-partition index at partition index i

    Base.append!(o, cp[ni][m])   # complete partition by appending it to a

    return o

end

function _restricted_partitions(o::Vector{Int}, n::Int, np::Int, cp::Vector{Vector{Vector{Int}}})

    oo = [o]

    for p = 1:np-1
        i = Base.findlast(x -> x > 1, oo[p])
        Base.append!(oo, [_partition(oo[p], n, i, cp)])
    end

    return oo

end

@doc raw"""
    integer_partitions(n [[[,m]; transpose=false], count=false])

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



