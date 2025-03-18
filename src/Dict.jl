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
#                            Dict.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                         dictSuperscript
# ------------------------------------------------------------------------------

@doc raw"""
    dictSuperscript 

```
julia> dictSuperscript
Dict{Char, Char} with 38 entries:
  'n' => 'ⁿ'
  'f' => 'ᶠ'
  'w' => 'ʷ'
  '1' => '¹'
  'd' => 'ᵈ'
  'e' => 'ᵉ'
  ⋮   => ⋮
```
#### Example:
```
julia> get(dictSuperscript, 'n', "not found")
'ⁿ': Unicode U+207F (category Lm: Letter, modifier)
```
"""
dictSuperscript = Dict(
    
    '-' => '⁻',
    '0' => '⁰',
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
    '/' => 'ᐟ', 
    
    'a' => 'ᵃ',
    'b' => 'ᵇ', 
    'c' => 'ᶜ',
    'd' => 'ᵈ',
    'e' => 'ᵉ',
    'f' => 'ᶠ',
    'g' => 'ᵍ',
    'h' => 'ʰ',
    'i' => 'ⁱ',
    'j' => 'ʲ',
    'k' => 'ᵏ',
    'l' => 'ˡ',
    'm' => 'ᵐ',
    'n' => 'ⁿ',
    'o' => 'ᵒ',
    'p' => 'ᵖ',
    'q' => '𐞥',
    'r' => 'ʳ',
    's' => 'ˢ',
    't' => 'ᵗ',
    'u' => 'ᵘ',
    'v' => 'ᵛ',
    'w' => 'ʷ',
    'x' => 'ˣ',
    'y' => 'ʸ',
    'z' => 'ᶻ'
    
    )

# ------------------------------------------------------------------------------
#                         dictSubscript
# ------------------------------------------------------------------------------

@doc raw"""
    dictSuperscript 

```
julia> dictSubscript
Dict{Char, Char} with 25 entries:
  'n' => 'ₙ'
  '1' => '₁'
  'e' => 'ₑ'
  '7' => '₇'
  '6' => '₆'
  'o' => 'ₒ'
  ⋮   => ⋮
```
#### Example:
```
julia> get(dictSubscript, 'n', "not found")
'ₙ': Unicode U+2099 (category Lm: Letter, modifier)
```
"""
dictSubscript = Dict(
    
    '-' => Char(0x208B),
    '0' => Char(0x2080),
    '1' => Char(0x2081),
    '2' => Char(0x2082),
    '3' => Char(0x2083),
    '4' => Char(0x2084),
    '5' => Char(0x2085),
    '6' => Char(0x2086),
    '7' => Char(0x2087),
    '8' => Char(0x2088),
    '9' => Char(0x2089),
    '/' => Char('⸝'),
    
    'a' => Char(0x2090),
    'e' => Char(0x2091),
    'h' => Char(0x2095),
    'k' => Char(0x2096),
    'l' => Char(0x2097),
    'm' => Char(0x2098),
    'n' => Char(0x2099),
    'o' => Char(0x2092),
    'p' => Char(0x209A),
    'r' => Char(0x1D63),
    's' => Char(0x209B),
    't' => Char(0x209C),
    'x' => Char(0x2093)
    
    )


# ------------------------------------------------------------------------------
#                         dictUndoSmall
# ------------------------------------------------------------------------------

@doc raw"""
    dictUndoSmall

```
julia> dictUndoSmall
Dict{Char, Char} with 62 entries:
  'ᵉ' => 'e'
  'ₓ' => 'x'
  'ₖ' => 'k'
  'ⁿ' => 'n'
  '𐞥' => 'q'
  'ₕ' => 'h'
  ⋮   => ⋮
```
#### Example:
```
julia> get(dictUndoSmall, 'ⁿ', "not found")
'n': ASCII/Unicode U+006E (category Ll: Letter, lowercase)
```
"""

dictUndoSmall = Dict(

    '₋' => '-', 
    '₀' => '0',
    '₁' => '1',
    '₂' => '2',
    '₃' => '3',
    '₄' => '4',
    '₅' => '5',
    '₆' => '6',
    '₇' => '7',
    '₈' => '8',
    '₉' => '9',
    '⸝' => '/',

    'ₐ' => 'a',
    'ₑ' => 'e',
    'ₕ' => 'h',
    'ₖ' => 'k',
    'ₗ' => 'l',
    'ₘ' => 'm',
    'ₙ' => 'n',
    'ₒ' => 'o',
    'ₚ' => 'p',
    'ᵣ' => 'r',
    'ₛ' => 's',
    'ₜ' => 't',
    'ₓ' => 'x',    

    '⁻' => '-',
    '⁰' => '0',
    '⁹' => '1',
    '²' => '2',
    '³' => '3',
    '⁴' => '4',
    '⁵' => '5',
    '⁶' => '6',
    '⁷' => '7',
    '⁸' => '8',
    '⁹' => '9',  
    'ᐟ' => '/',
 
    'ᵃ' => 'a',
    'ᵇ' => 'b', 
    'ᶜ' => 'c',
    'ᵈ' => 'd',
    'ᵉ' => 'e',
    'ᶠ' => 'f',
    'ᵍ' => 'g',
    'ʰ' => 'h',
    'ⁱ' => 'i',
    'ʲ' => 'j',
    'ᵏ' => 'k',
    'ˡ' => 'l',
    'ᵐ' => 'm',
    'ⁿ' => 'n',
    'ᵒ' => 'o',
    'ᵖ' => 'p',
    '𐞥' => 'q',
    'ʳ' => 'r',
    'ˢ' => 's',
    'ᵗ' => 't',
    'ᵘ' => 'u',
    'ᵛ' => 'v',
    'ʷ' => 'w',
    'ˣ' => 'x',
    'ʸ' => 'y',
    'ᶻ' => 'z'
    
    )