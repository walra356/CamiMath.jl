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

# ==============================================================================
#                            Pochhammer.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#               pochhammer(x::T, p::Int) where T<:Real
# ------------------------------------------------------------------------------

@doc raw"""
    pochhammer(x::T, p::Int) where T<:Real

Pochhammer symbol ``(x)_{p}`` for *non-negative* integer `p`,
```math
(x)_{p}=\begin{cases}
1 & p=0\\
x(x+1)(x+2)⋯(x+p-1) & p>0
\end{cases}
```
Note that ``(x)_{p}=0`` for ``x=0,-1,⋯\ -(p-1)``
#### Examples:
```
julia> x = [-4,-3,-2,-1, 0, 1, 2 , 3, 4];

julia> pochhammer.(x, 5) == [0, 0, 0, 0, 0, 120, 720, 2520, 6720]
true

julia> pochhammer.(x, 0) == [1, 1, 1, 1, 1, 1, 1, 1, 1]
true

julia> o = [pochhammer.([x for x=0:-1:-p],p) for p=0:5]
6-element Vector{Vector{Int64}}:
 [1]
 [0, -1]
 [0, 0, 2]
 [0, 0, 0, -6]        
 [0, 0, 0, 0, 24]     
 [0, 0, 0, 0, 0, -120]

julia>  o = [pochhammer.([x for x=0:p],p) for p=0:5]
6-element Vector{Vector{Int64}}:
 [1]
 [0, 1]
 [0, 2, 6]
 [0, 6, 24, 60]
 [0, 24, 120, 360, 840]
 [0, 120, 720, 2520, 6720, 15120]

julia> x = Rational{BigInt}(-1, 50);

julia> pochhammer(x, 20)
-21605762356630090481082546653745369902321614221999//9536743164062500000000000000000000
```
"""
function pochhammer(x::T, p::Int) where {T<:Real}

  p > 0 || return 1

  o = x

  for n = 1:p-1
    o *= (x + n)
  end

  return o

end
