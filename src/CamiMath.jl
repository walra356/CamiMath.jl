# SPDX-License-Identifier: MIT

# =============================================================================
#                               CamiMath.jl
#                              Jook Walraven
# =============================================================================

module CamiMath

export fwd
export bwd
export reg
export rev
export Type_IOP

export isforward
export isbackward
export isregular
export isreversed

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
export log10_characteristic
export log10_mantissa
export normalize_rationals
export numerators
export pascal_triangle
export pascal_next
export permutations_unique_count

export laguerreL
export laguerre_polynom
export generalized_laguerre_polynom
export generalized_laguerreL
export lagrange_polynom

export pochhammer
export polynomial
export polynom_power
export polynom_product
export polynom_product_expansion
export texp

export istriangle
export triangle_coefficient
export threeJsymbol
export CGC

include("Singleton.jl")
include("Bernoulli.jl")
include("Divisor.jl")
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
include("Lagrange.jl")
include("Toolbox.jl")
include("VectorCoupling.jl")

end
