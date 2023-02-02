# SPDX-License-Identifier: MIT

# ==============================================================================
#                        Bigfactorial.jl
# ==============================================================================

@doc raw"""
    bigfactorial(n::Int [; msg=true])

The product of all *positive* integers less than or equal to `n`,
```math
n!=n(n-1)(n-2)⋯1.
```
By definition
```math
0!=1
```
For *negative* integers the factorial is zero. Integer-overflow protection: 
for `n` > 20 the output is autoconverted to `BigInt`. 
By default the capture message is activated: 
"Warning (IOP): bigfactorial converted to BigInt".
#### Examples:
```
julia> bigfactorial(20) == factorial(20)
true

julia> bigfactorial(21)
Warning (IOP): bigfactorial converted to BigInt
51090942171709440000

julia> bigfactorial(21; msg=false)
51090942171709440000

julia> factorial(21)
Warning (IOP): bigfactorial converted to BigInt
```
"""
function bigfactorial(n::T; msg=true) where {T<:Integer}

    str = "Warning (IOP): bigfactorial converted to BigInt"

    T == Int || return factorial(n)

    n = Int(n)
    nc = 20

    o = n > nc ? factorial(big(n)) : o = factorial(n)

    msg && n > nc && println(str)

    return o

end