# SPDX-License-Identifier: MIT

# ==============================================================================
#                             Toolbox.jl
#                      Jook Walraven - 23-2-2023
# ==============================================================================

# ------------------------------------------------------------------------------
#         Type_IOP(n::Integer, nc::Integer [, a [; fnam="" [; msg=true]]])
# ------------------------------------------------------------------------------

@doc raw"""
    Type_IOP(n::Integer, nc::Integer [, a [; fnam="" [; msg=true]]])

`BigInt` if `n` is a `BigInt` or `n > nc`, otherwise `Int`; `a` is an 
auxiliary second variable.

- `fnam` : function name

- `msg` : integer-overflow protection (IOP) - warning on activation 
#### Examples:
```
julia> Type_IOP(1, 1)
Int64

julia> Type_IOP(big(1), 1)
BigInt

julia> Type_IOP(2, 1)
BigInt

julia> Type_IOP(1, 1; fnam="test")
Int64

julia> Type_IOP(2, 1, 0; fnam="test")
 IOP capture at test(2, 0): output converted to BigInt
BigInt
```
"""
function Type_IOP(n::Integer, nc::Integer, a=nothing; fnam="", msg=true)

    warning = " output converted to BigInt\n"

    if n isa BigInt
        return BigInt
    else
        isempty(fnam) ? nothing :
        n ≤ nc ? nothing :
        !msg ? nothing :
        isnothing(a) ? print(" IOP capture at " * fnam * "($n):" * warning) :
        print(" IOP capture at " * fnam * "($n, $a):" * warning)

        return n ≤ nc ? Int : BigInt
    end

end