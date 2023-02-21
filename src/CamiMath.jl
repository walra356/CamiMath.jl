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
export faulhaber_polynomial1
export faulhaber_summation
export fibonacci
export harmonicNumber
export istriangle
export canonical_partitions
export integer_partitions
export log10_characteristic_power
export log10_mantissa
export normalize_rationals
export numerators
export pascal_triangle
export pascal_next
export permutations_unique_count

export laguerre_polynoms
export laguerre_polynom
export generalized_laguerre_polynoms

export pochhammer
export polynomial
export polynom_power
export polynom_product
export polynom_product1
export polynom_product_expansion
export texp
export triangle_coefficient

include("Bernoulli.jl")
include("Divisor.jl")
include("Exponential.jl")
include("Factorial.jl")
include("Faulhaber.jl")
include("Fibonacci.jl")
include("HarmonicNumber.jl")
include("Laguerre.jl")
include("Partition.jl")
include("Permutations.jl")
include("Pascal.jl")
include("Pochhammer.jl")
include("Polynomials.jl")
include("Triangle.jl")

end
