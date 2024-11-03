# SPDX-License-Identifier: MIT

# author: Jook Walraven - 3-11-2024

# ==============================================================================
#                              Structs.jl
# ==============================================================================

@doc raw"""
    fwd

Singleton type indicating ``forward``
"""
struct fwd
end

@doc raw"""
    bwd

Singleton type indicating ``backward``
"""
struct bwd
end

@doc raw"""
    reg

Singleton type indicating ``regular``
"""
struct reg
end

@doc raw"""
    rev

Singleton type indicating ``reverse``
"""
struct rev
end

# ============================= isforward(notation) ===========================

@doc raw"""
function isforward(sense)

Boolean status of `sense`, with options: `fwd` (forward) and `bwd` (backward).
#### Example:
```
julia> isforward(fwd)
true
```
"""
function isforward(sense)

strErr = "Error: invalid sense (options: fwd, bwd)"

return sense === fwd ? true : sense === bwd ? false : error(strErr)

end

# ============================= isregular(ordering) ============================
@doc raw"""
function isregular(sense::Type)

Boolean status of `sense`, with options: `reg` (regular) and `rev` (reversed).
#### Example:
```
julia> isregular(reg)
true
```
"""
function isregular(sense)

strErr = "Error: invalid sense (options: reg, rev)" 

return sense === reg ? true : sense === rev ? false : error(strErr)

end
# ============================= End ===========================
