# SPDX-License-Identifier: MIT

# ==============================================================================
#                            Permutations.jl
#                       Jook Walraven - 17-2-2023
# ==============================================================================

# ------------------------------------------------------------------------------
#         permutations_unique_count(p::Vector{Vector{Int}}, i::Int)
# ------------------------------------------------------------------------------

@doc raw"""
    permutations_unique_count(p::Vector{Vector{Integer}}, i::Int)

Number of unique permutations of the subarray ``p[i]``.
#### Example:
```
p = [[1,2,3],[2,3,1,4,3]]
permutations_unique_count(p,2)
 60
```
"""
function permutations_unique_count(p::Vector{Vector{Int}}, i::Int)

    o = CamiMath.bigfactorial(length(p[i]))
    d = Dict([(n, Base.count(x -> x == n, p[i])) for n ∈ Base.unique(p[i])])

    for j ∈ Base.eachindex(Base.unique(p[i]))
        o = o ÷ CamiMath.bigfactorial(d[Base.unique(p[i])[j]])
    end

    return o

end