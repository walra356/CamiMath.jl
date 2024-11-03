# SPDX-License-Identifier: MIT

# author: Jook Walraven - 3-11-2024

# ==============================================================================
#                              Lagrange.jl
# ==============================================================================

@doc raw"""
    lagrange_polynom(f::Vector{T}, start::Int, stop::Int [, sense=fwd]) where T<:Real
    lagrange_polynom(f::Vector{T}, itr::UnitRange [, sense=fwd]) where T<:Real

The coefficients of the polynomial of degree ``d = stop-start`` running through
``d+1`` subsequent points, ``f[n::n+d]``, of the tabulated regular function ``f[n]``. 
The corresponding polynomial expansion is most accurate near ``start``/``stop`` (for ``fwd/bwd``). 
#### Examples:
```
julia> ftab = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];

julia> coords_fwd = lagrange_polynom(ftab, 1, 4, fwd); println(coords_fwd)
[0.0, 0.0, 1.0, 0.0]

julia> coords_bwd = lagrange_polynom(ftab, 1, 4, bwd); println(coords_bwd)
[9.0, 6.0, 1.0, -4.440892098500626e-16]

julia> f(x) = polynomial(coords_fwd, x); 

julia> f(0.0), f(0.5), f(1.0), f(2.0), f(3.0)
(0.0, 0.25, 1.0, 4.0, 9.0)

julia> f(x) = polynomial(coords_bwd, x); 

julia> f(-3.0), f(-2.5), f(-2.0), f(-1.0), f(0.0)
(1.199040866595169e-14, 0.25000000000000694, 1.0000000000000036, 4.0, 9.0)
```
"""
function lagrange_polynom(f::Vector{T}, start::Int, stop::Int, sense=fwd) where T<:Real

    k = stop-start

    coords = isforward(sense) ? inv([T(m)^j for m=0:k, j=0:k])*f[start:stop] : 
                                  inv([T(-m)^j for m=0:k, j=0:k])*f[stop:-1:start]
        
    return coords
    
end
function lagrange_polynom(f::Vector{T}, itr::UnitRange, sense=fwd) where T<:Real
        
    return lagrange_polynom(f, itr.start, itr.stop, sense)
    
end