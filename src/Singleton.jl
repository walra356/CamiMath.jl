# SPDX-License-Identifier: MIT

# author: Jook Walraven - 3-11-2024

# ==============================================================================
#                              Singleton.jl
# ==============================================================================

@doc raw"""
    fwd

Singleton type indicating `forward` sense
"""
struct fwd
end

@doc raw"""
    bwd

Singleton type indicating `backward` sense
"""
struct bwd
end

@doc raw"""
    reg

Singleton type indicating `regular` ordering
"""
struct reg
end

@doc raw"""
    rev

Singleton type indicating `reverse` ordering
"""
struct rev
end