# SPDX-License-Identifier: MIT

# Copyright (c) 2023 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# =============================================================================
#                               CamiMath.jl
# =============================================================================

module CamiMath

import Documenter
#import DocumenterInterLinks

export sup
export sub
export frac
export strRational
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
export polynom
export polynom_power
export polynom_product
export polynom_product_expansion
export texp

export istriangle
export triangle_coefficient
export threeJsymbol
export CGC

include("Polynomials.jl")
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
include("Lagrange.jl")
include("Toolbox.jl")
include("VectorCoupling.jl")

end
