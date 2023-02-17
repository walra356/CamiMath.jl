# SPDX-License-Identifier: MIT

# ==============================================================================
#                            Pochhammer.jl
#                      Jook Walraven - 10-2-2023
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
