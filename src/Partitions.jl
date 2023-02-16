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



