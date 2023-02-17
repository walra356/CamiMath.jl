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
export canonical_partitions
export integer_partitions
export normalize_rationals
export numerators
export pascal_triangle
export pascal_next
export pochhammer

include("Bernoulli.jl")
include("Divisor.jl")
include("Factorial.jl")
include("Faulhaber.jl")
include("Fibonacci.jl")
include("HarmonicNumber.jl")
include("Partition.jl")
Include("Pochhammer.jl")
include("Pascal.jl")

end
