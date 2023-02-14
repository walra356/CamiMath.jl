# SPDX-License-Identifier: MIT

# =============================================================================
#                               CamiMath.jl
# =============================================================================

module CamiMath

export bernoulliB
export bigfactorial
export faulhaber_polynom
export faulhaber_polynomial
export faulhaber_summation
export harmonicNumber
export pascal_triangle
export pascal_next_row

include("Bernoulli.jl")
include("Factorial.jl")
include("Faulhaber.jl")
include("HarmonicNumber.jl")
include("Pascal.jl")

end
