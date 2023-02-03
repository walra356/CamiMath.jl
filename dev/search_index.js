var documenterSearchIndex = {"docs":
[{"location":"#CamMath.jl","page":"Home","title":"CamMath.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Mathematics library with integer-overload protection","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Bernoulli-number","page":"Home","title":"Bernoulli number","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"bernoulliB(n::T; msg=true) where {T<:Integer}   \nbernoulliB_array(nmax::T; msg=true) where {T<:Integer}","category":"page"},{"location":"#CamMath.bernoulliB-Tuple{T} where T<:Integer","page":"Home","title":"CamMath.bernoulliB","text":"bernoulliB(n::T [; msg=true]) where {T<:Integer}\n\nBernoulli numbers of index n are defined by the recurrence relation\n\n    B_n = - frac1n+1sum_k=0^n-1frac(n+1)k(n+1-k)B_k\n\nwith B_0=1 and B_1=-12. Starting at B_0 is called the even index  convention (B_2n+1=0 rmfor n1).\n\nInteger-overflow protection (IOP): for n > 35 the output is autoconverted to  Rational{BigInt}. By default the capture message is activated:  \"Warning: bernoulliB converted to Rational{BigInt}\".\n\nExamples:\n\njulia> o = [bernoulliB(n) for n=0:5]; println(o)\nRational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]\n\njulia> bernoulliB(60)\nWarning: bernoulliB converted to Rational{BigInt}\n-1215233140483755572040304994079820246041491//56786730\n\njulia> n = 60;\njulia> bernoulliB(n) == bernoulliB_array(n)[end]             \ntrue\n\n\n\n\n\n","category":"method"},{"location":"#CamMath.bernoulliB_array-Tuple{T} where T<:Integer","page":"Home","title":"CamMath.bernoulliB_array","text":"bernoulliB_array(nmax::T [; msg=true]) where {T<:Integer}\n\nBernoulli number array B_0cdots B_nmax, where nmax is the index of  the highest Bernoulli number of the array (NB.: not the array length).\n\nExamples:\n\njulia> o = bernoulliB_array(8); println(o)\nRational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30]\n\njulia> o = bernoulliB_array(big(8)); println(o)\nRational{BigInt}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30]\n\njulia> n = 60; msg = false;\njulia>  bernoulliB(n; msg) == bernoulliB_array(n; msg)[end]            \ntrue\n\n\n\n\n\n","category":"method"},{"location":"#Factorial","page":"Home","title":"Factorial","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"bigfactorial(n::T; msg=true) where {T<:Integer}","category":"page"},{"location":"#CamMath.bigfactorial-Tuple{T} where T<:Integer","page":"Home","title":"CamMath.bigfactorial","text":"bigfactorial(n::Int [; msg=true])\n\nThe product of all positive integers less than or equal to n,\n\nn=n(n-1)(n-2)1\n\nIn addition 0=1 by definition. For negative integers the factorial is zero. Integer-overflow protection  (IOP): for n > 20 the output is autoconverted to BigInt.\n\nExamples:\n\njulia> bigfactorial(20) == factorial(20)\ntrue\n\njulia> bigfactorial(21)\nWarning: bigfactorial converted to BigInt\n51090942171709440000\n\njulia> bigfactorial(21; msg=false)\n51090942171709440000\n\njulia> factorial(21)\nERROR: OverflowError: 21 is too large to look up in the table; consider using \n`factorial(big(21))` instead\n\n\n\n\n\n","category":"method"},{"location":"#Faulhaber-polynomial","page":"Home","title":"Faulhaber polynomial","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"faulhaber_polynom(k::Int; T=Int)\nfaulhaber_polynomial(n::T, p::T) where {T<:Integer}\nfaulhaber_summation(n::Int, p::Int)","category":"page"},{"location":"#CamMath.faulhaber_polynom-Tuple{Int64}","page":"Home","title":"CamMath.faulhaber_polynom","text":"faulhaber_polynom(p::T) where {T<:Integer}\n\nVector representation of the coefficients of the faulhaber_polynomial  of degree p \n\n   c=c_0 c_p\n\nwith vector elements\n\n    c_0=0  rmand c_j=frac1pbinom pp-jB_p-j\n\nwhere j 1 p. The B_0 B_p-1 are Bernoulli numbers (but with B_1=+frac12 rather than -frac12).\n\nExample:\n\nfaulhaber_polynom(6)\n7-element Vector{Rational{Int64}}:\n  0//1\n  0//1\n -1//12\n  0//1\n  5//12\n  1//2\n  1//6\n\n\n\n\n\n","category":"method"},{"location":"#CamMath.faulhaber_polynomial-Union{Tuple{T}, Tuple{T, T}} where T<:Integer","page":"Home","title":"CamMath.faulhaber_polynomial","text":"faulhaber_polynomial(n::T, p::T) where {T<:Integer}\n\nFaulhaber polynomial of degree p \n\n    F(np)=sum_j=0^pc_jn^j\n\nwhere the coefficients are contained in the coefficient vector  faulhaber_polynom.\n\nExample:\n\n\n\n\n\n\n\n","category":"method"},{"location":"#CamMath.faulhaber_summation-Tuple{Int64, Int64}","page":"Home","title":"CamMath.faulhaber_summation","text":"faulhaber_summation(n::T, p::T) where {T<:Integer}\n\nSum of powers of natural numbers 1 n as given by the Faulhaber formula\n\n    sum_k=1^nk^p=F(np+1)\n\nwhere F(np) is the faulhamer_polynomial of degree p.\n\nExamples:\n\nfaulhaber_summation(5,1)\n 15\n\nfaulhaber_summation(3,60; T=BigInt)\n  42391158276369125018901280178\n\n\n\n\n\n","category":"method"},{"location":"#HarmonicNumber","page":"Home","title":"HarmonicNumber","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"harmonicNumber(n::T; msg=true) where {T<:Integer}\nharmonicNumber_array(nmax::T; msg=true) where {T<:Integer}\nharmonicNumber(n::T, p::Int; msg=true) where {T<:Integer}","category":"page"},{"location":"#CamMath.harmonicNumber-Tuple{T} where T<:Integer","page":"Home","title":"CamMath.harmonicNumber","text":"harmonicNumber(n::T [; msg=true]) where {T<:Integer}\n\nSum of the reciprocals of the first n natural numbers\n\n    H_n=sum_k=1^nfrac1k\n\nInteger-overflow protection: for n > 46 the output is autoconverted to Rational{BigInt}. By default the capture message is activated:  \"Warning: harmonicNumber autoconverted to Rational{BigInt}\". \n\nExamples:\n\njulia> o = harmonicNumber_array(9); println(o)\nRational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280, 7129//2520]\n\njulia> o = [harmonicNumber(46; msg=true)]; println(o)\nRational{Int64}[5943339269060627227//1345655451257488800]\n\njulia> o = [harmonicNumber(47; msg=true)]; println(o)\nWarning: harmonicNumber autoconverted to Rational{BigInt}\nRational{BigInt}[282057509927739620069//63245806209101973600]\n\njulia> harmonicNumber(12) == harmonicNumber(12, 1)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#CamMath.harmonicNumber_array-Tuple{T} where T<:Integer","page":"Home","title":"CamMath.harmonicNumber_array","text":"harmonicNumber_array(nmax::T [; msg=true]) where {T<:Integer}\n\nSum of the reciprocals of the first n natural numbers\n\n    H_n=sum_k=1^nfrac1k\n\nInteger-overflow protection: for n > 46 the output is autoconverted to Rational{BigInt}. By default the capture message is activated:  \"Warning: harmonicNumber autoconverted to Rational{BigInt}\". \n\nExamples:\n\njulia> o = harmonicNumber_array(9); println(o)\nRational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280, 7129//2520]\n\njulia> o = [harmonicNumber(46; msg=true)]; println(o)\nRational{Int64}[5943339269060627227//1345655451257488800]\n\njulia> o = [harmonicNumber(47; msg=true)]; println(o)\nWarning: harmonicNumber autoconverted to Rational{BigInt}\nRational{BigInt}[282057509927739620069//63245806209101973600]\n\njulia> harmonicNumber(12) == harmonicNumber(12, 1)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#CamMath.harmonicNumber-Union{Tuple{T}, Tuple{T, Int64}} where T<:Integer","page":"Home","title":"CamMath.harmonicNumber","text":"harmonicNumber(n::T, p::Int [; msg=true]) where {T<:Integer}\n\nSum of the p^th power of reciprocals of the first n numbers\n\n    H_np=sum_k=1^nfrac1k^p\n\nInteger-overflow protection: the output is autoconverted to Rational{BigInt} when required. By default the capture message is activated:  \"Warning: harmonicNumber autoconverted to Rational{BigInt}\". \n\nExamples:\n\njulia> o = [harmonicNumber(46,1; msg=true)]; println(o)\nRational{Int64}[5943339269060627227//1345655451257488800]\n\njulia> o = [harmonicNumber(47,1; msg=true)]; println(o)\nWarning: harmonicNumber autoconverted to Rational{BigInt}\"\nRational{BigInt}[280682601097106968469//63245806209101973600]\n\njulia> o = [harmonicNumber(47,1)]; println(o)\nRational{BigInt}[280682601097106968469//63245806209101973600]\n\nharmonicNumber(12, -3) == faulhaber_summation(12, 3)\n  true\n\n\n\n\n\n","category":"method"},{"location":"#Pascal-triangle","page":"Home","title":"Pascal triangle","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pascal_triangle(nmax::T) where {T<:Integer}\npascal_next(a::Vector{Int})","category":"page"},{"location":"#CamMath.pascal_triangle-Tuple{T} where T<:Integer","page":"Home","title":"CamMath.pascal_triangle","text":"pascal_triangle(nmax [, T=Int])\n\nPascal triangle of binomial coefficients binomnk for n=0 1 nmax\n\nExample:\n\npascal_triangle(5)\n6-element Vector{Vector{Int64}}:\n [1]\n [1, 1]\n [1, 2, 1]\n [1, 3, 3, 1]\n [1, 4, 6, 4, 1]\n [1, 5, 10, 10, 5, 1]\n\n\n\n\n\n","category":"method"},{"location":"#CamMath.pascal_next-Tuple{Vector{Int64}}","page":"Home","title":"CamMath.pascal_next","text":"pascal_next(nmax)\n\nNext row of Pascal triangle\n\nExample:\n\na = [1, 4, 6, 4, 1]\npascal_next(a)\n [1, 5, 10, 10, 5, 1]\n\n\n\n\n\n","category":"method"}]
}
