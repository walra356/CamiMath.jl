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
export harmonicNumber
export harmonicNumber_array

include("Bernoulli.jl")
include("Factorial.jl")
include("Faulhaber.jl")
include("HarmonicNumber")

end
