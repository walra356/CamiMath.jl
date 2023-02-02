# SPDX-License-Identifier: MIT

# =============================================================================
#                               CamMath.jl
# =============================================================================

module CamMath

export bernoulliB
export bernoulliB_array
export bigfactorial
export faulhaber_polynom
export faulhaber_summation
export harmonicNumber
export harmonicNumber_array
export pascal_triangle
export pascal_next

include("Bernoulli.jl")
include("Factorial.jl")
include("Faulhaber.jl")
include("HarmonicNumber.jl")
include("Pascal.jl")

end
