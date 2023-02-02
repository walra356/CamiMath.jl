# SPDX-License-Identifier: MIT

# =============================================================================
#                               CamMath.jl
# =============================================================================

module CamMath

export bernoulliB
export bernoulliB_array
export bigfactorial

include("BernoulliB.jl")
include("Bigfactorial.jl")

end
