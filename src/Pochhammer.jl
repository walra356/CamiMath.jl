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

Pochhammer symbol ``(x)_{p}`` for integral ``p``,
```math
(x)_{p}=\begin{cases}
1 & p=0\\
x(x+1)(x+2)⋯(x+p-1) & p>0
\end{cases}
```
Note that ``(x)_{p}=0`` for ``x=0,-1,⋯\ -(p-1)``
#### Examples:
```
x = [-4,-3,-2,-1, 0, 1, 2 , 3, 4]
pochhammer.(x,5) == [0, 0, 0, 0, 0, 120, 720, 2520, 6720]
  true
pochhammer.(x,0) == [1, 1, 1, 1, 1, 1, 1, 1, 1]
  true
o = [pochhammer.([x for x=0:-1:-p],p) for p=0:5]
println("non-positive integer x = 0,⋯\ -p:")
for p=0:5
    println("p = $p: $(o[p+1])")
end
  non-positive integer x = 0,⋯\ -p:
  p = 0: [1]
  p = 1: [0, -1]
  p = 2: [0, 0, 2]
  p = 3: [0, 0, 0, -6]
  p = 4: [0, 0, 0, 0, 24]
  p = 5: [0, 0, 0, 0, 0, -120]
 o = [pochhammer.([x for x=0:p],p) for p=0:5]
 println("non-negative integer x = 0,⋯\  p:")
 for p=0:5
     println("p = $p: $(o[p+1])")
 end
   non-negative integer x = 0,⋯\  p:
   p = 0: [1]
   p = 1: [0, 1]
   p = 2: [0, 2, 6]
   p = 3: [0, 6, 24, 60]
   p = 4: [0, 24, 120, 360, 840]
   p = 5: [0, 120, 720, 2520, 6720, 15120]
x = -1//50
pochhammer(x,20)
  OverflowError: -1491212300990613201 * 449 overflowed for type Int64
x = convert(Rational{BigInt}, -1//50)
pochhammer(x,20)
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
