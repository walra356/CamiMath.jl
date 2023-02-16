
# SPDX-License-Identifier: MIT

# author: Jook Walraven - 15-2-2023

# ==============================================================================
#                               Fibonacci.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#            fibonacci(n::Integer [[; arr=false], msg=true])
# ------------------------------------------------------------------------------

function _fibonacci_next!(o, k)

    while k > 0
        Base.append!(o, o[end-1] + o[end])
        k -= 1
    end

    return o

end

# ..............................................................................
@doc raw"""
    fibonacci(n::Integer [[; arr=false], msg=true])

The sequence of integers,  ``F_0,⋯\ F_{nmax}``, in which each element is 
the sum of the two preceding ones, 
```math
    F_n = F_{n-1}+F_{n-2}.
```
with ``F_1=1`` and ``F_0=0``. 
Integer-overflow protection: for `n > 92` the output is autoconverted to BigInt. 
By default the capture message is activated: "Warning: fibonacciF autoconverted to BigInt". 
#### Example:
```
julia> fibonacciF_array(20)
21-element Vector{Int64}:
   0
   1
   1
   2
   3
   5
   ⋮
1597
2584
4181
julia> fibonacciF(100)
Warning: fibonacciF autoconverted to BigInt
354224848179261915075
```
"""
function fibonacci(n::Integer; arr=false, msg=true)

    o = (
        1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
        1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393,
        196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887,
        9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141,
        267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073,
        4807526976, 7778742049, 12586269025, 20365011074, 32951280099,
        53316291173, 86267571272, 139583862445, 225851433717, 365435296162,
        591286729879, 956722026041, 1548008755920, 2504730781961,
        4052739537881, 6557470319842, 10610209857723, 17167680177565,
        27777890035288, 44945570212853, 72723460248141, 117669030460994,
        190392490709135, 308061521170129, 498454011879264, 806515533049393,
        1304969544928657, 2111485077978050, 3416454622906707, 5527939700884757,
        8944394323791464, 14472334024676221, 23416728348467685,
        37889062373143906, 61305790721611591, 99194853094755497,
        160500643816367088, 259695496911122585, 420196140727489673,
        679891637638612258, 1100087778366101931, 1779979416004714189,
        2880067194370816120, 4660046610375530309, 7540113804746346429
    )

    nc = 92

    T = n > nc ? BigInt : typeof(n)

    n = convert(Int, n)

    if n < 0
        throw(DomainError(n))
    elseif n == 0
        return [T(0)]
    elseif n > nc
        str = "IOP capture: fibonaci($n) converted to BigInt"
        msg && T == Int && println(str)
    end

    if arr # -------------------------------------------------------------------
        if n ≤ nc
            return T[o[i] for i = 1:n]
        else
            o = collect(BigInt, o)
            return _fibonacci_next!(o, n - nc)
        end
    else # ---------------------------------------------------------------------
        if n ≤ nc
            return T.(o)[n]
        else
            o = collect(BigInt, (o[end-1], o[end]))
            return _fibonacci_next!(o, n - nc)[end]
        end
    end # ----------------------------------------------------------------------

end