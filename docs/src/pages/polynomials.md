# Polynomials

A set of polynomials can be regarded as the elements of a vector space. As an example 
we consider the set of all *real* polynomials of degree ``d``,
```math
P_α(x) = α_0 + α_1 x + ⋯ + α_d x^d.
```
These are maps ``P_α:\mathbb{\mathbb{R\rightarrow\mathbb{R}}}`` that 
satisfy the group operation 'addition of polynomials' because the sum of two 
polynomials of degree ``d`` is again a polynomial of degree ``d``,
```math
(P_α + P_β)(x) ≡ (α_0 + β_{0})+(α_1 + β_1)x + ⋯ + (α_d + β_d) x^d = P_α(x) + P_β(x),
```
and remains a polynomial of degree ``d`` under 'scalar multiplication',
```math
(λ P_α)(x) ≡ P_{\lambdaα}(x)=\lambdaα_0+\lambdaα_1 x + ⋯ + \lambdaα_d x^d = λ P_α(x).
```
The 'zero element' of the vector space is the polynomial ``P_α(x)=0`` and the 'inverse 
element' of the element ``P_α(x)`` is the polynomial ``-P_α(x)``. Also the 'associative' 
and 'distributive' properties are easily verified. 

Hence, the set of all *real* polynomials of degree ``d`` defines  a vector space (of 
dimension ``d + 1``) over the field ``\mathbb{R}`` and the polynomials 
``1,x,x^{2},\cdots x^d`` represent a *basis*. The coefficients ``α_0, ⋯, α_d``,
represent the *coordinates* of the vector ``P_α`` with respect to this basis.

In `CamiMath` we define a polynomial by specifying the coordinate vector [`polynom`](@ref).

### polynom

```@docs
polynom
```

### polynom power

```@docs
polynom_power(polynom, power::Int)
```

### polynom product

```@docs
polynom_product(polynom1, polynom2)
```

### polynom product expansion

```@docs
polynom_product_expansion(polynom1, polynom2, p::Int)
```

### polynomial

```@docs
polynomial(polynom, x::T; deriv=0) where T<:Real
```

