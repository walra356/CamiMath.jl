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

- `msg` : integer-overflow protection (IOP) - warning on activation  (for `n > 20`) 
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
function bigfactorial(n::Integer; msg=true)

    str = "IOP capture: bigfactorial($n) converted to BigInt"

    typeof(n) == Int || return factorial(n)

    n = Int(n)
    nc = 20

    o = n > nc ? factorial(big(n)) : factorial(n)

    msg && n > nc && println(str)

    return o

end