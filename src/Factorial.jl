# SPDX-License-Identifier: MIT

# ==============================================================================
#                              Factorial.jl
# ==============================================================================

@doc raw"""
    bigfactorial(n::Int [; msg=true])

The product of all *positive* integers less than or equal to `n`,
```math
n!=n(n-1)(n-2)⋯1.
```
In addition ``0!=1`` by definition.
For *negative* integers the factorial is zero. 
Integer overflow protection (IOP): on integer overflow the output is converted 
to Rational{BigInt}. By default the IOP capture message is activated.
#### Examples:
```
julia> bigfactorial(20) == factorial(20)
true

julia> bigfactorial(21)
IOP capture: bigfactorial(21) converted to BigInt
51090942171709440000

julia> bigfactorial(21; msg=false)
51090942171709440000

julia> factorial(21)
ERROR: OverflowError: 21 is too large to look up in the table; consider using 
`factorial(big(21))` instead
```
"""
function bigfactorial(n::T; msg=true) where {T<:Integer}

    str = "IOP capture: bigfactorial($n) converted to BigInt"

    T == Int || return factorial(n)

    n = Int(n)
    nc = 20

    o = n > nc ? factorial(big(n)) : o = factorial(n)

    msg && n > nc && println(str)

    return o

end