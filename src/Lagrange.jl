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
#                              Lagrange.jl
# ==============================================================================

@doc raw"""
    lagrange_polynom(f::Vector{T}, start::Int, stop::Int [, sense=fwd]) where T<:Real
    lagrange_polynom(f::Vector{T}, itr::UnitRange [, sense=fwd]) where T<:Real

The coefficients of the [`polynomial`](@ref) of degree ``d = ```stop`-`start` running 
through ``d+1`` subsequent points of the tabulated regular function ``f[n]``. 
For `sense` = [`fwd](@ref) these are the points ``f[n:n+d]``, for `sense` = [`bwd`](@ref) 
the points ``f[n-k:n]``.
The corresponding [`polynomial`](@ref) is most accurate near ``f[n]``. 
#### Examples:
```
julia> a = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];

julia> polynom = lagrange_polynom(a, 1, 4, fwd); println(polynom)
[0.0, 0.0, 1.0, 0.0]

julia> f(x) = polynomial(polynom, x); 

julia> f(0.0), f(0.5), f(1.0), f(2.0), f(3.0)
(0.0, 0.25, 1.0, 4.0, 9.0)

julia> polynom = lagrange_polynom(a, 1:4, bwd); println(polynom)
[9.0, 6.0, 1.0, -4.440892098500626e-16]

julia> f(x) = polynomial(polynom, x); 

julia> f(-3.0), f(-2.5), f(-2.0), f(-1.0), f(0.0)
(1.199040866595169e-14, 0.25000000000000694, 1.0000000000000036, 4.0, 9.0)
```
"""
function lagrange_polynom(f::Vector{T}, start::Int, stop::Int, sense=fwd) where T<:Real

    k = stop-start

    polynom = isforward(sense) ? inv([T(m)^j for m=0:k, j=0:k])*f[start:stop] : 
                                inv([T(-m)^j for m=0:k, j=0:k])*f[stop:-1:start]
        
    return polynom
    
end
function lagrange_polynom(f::Vector{T}, itr::UnitRange, sense=fwd) where T<:Real
        
    return lagrange_polynom(f, itr.start, itr.stop, sense)
    
end