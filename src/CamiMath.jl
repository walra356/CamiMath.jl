# SPDX-License-Identifier: MIT

# =============================================================================
#                               CamiMath.jl
#                              Jook Walraven
# =============================================================================

module CamiMath

export bernoulliB
export bigfactorial
export divisor
export faulhaber_polynom
export faulhaber_polynomial
export faulhaber_summation
export fibonacci
export harmonicNumber
export istriangle
export canonical_partitions
export integer_partitions
export normalize_rationals
export numerators
export pascal_triangle
export pascal_next
export permutations_unique_count
export pochhammer
export texp
export triangle_coefficient

include("Bernoulli.jl")
include("Divisor.jl")
include("Exponential.jl")
include("Factorial.jl")
include("Faulhaber.jl")
include("Fibonacci.jl")
include("HarmonicNumber.jl")
include("Partition.jl")
include("Permutations.jl")
include("Pochhammer.jl")
include("Pascal.jl")
include("Triangle.jl")

end
