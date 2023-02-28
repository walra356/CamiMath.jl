var documenterSearchIndex = {"docs":
[{"location":"man/library/#Bernoulli-number","page":"Library","title":"Bernoulli number","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"bernoulliB(n::Integer; arr=false, msg=true)","category":"page"},{"location":"man/library/#CamiMath.bernoulliB-Tuple{Integer}","page":"Library","title":"CamiMath.bernoulliB","text":"bernoulliB(n::Integer [; arr=false [, msg=true]])\n\nBernoulli numbers of index n are defined by the recurrence relation\n\n    B_n = - frac1n+1sum_k=0^n-1frac(n+1)k(n+1-k)B_k\n\nwith B_0=1 and B_1=-12. Including B_0 results in the even index  convention (B_2n+1=0 for n1).\n\narr : output in array format\nmsg : integer-overflow protection (IOP) - warning on activation \n\nExamples:\n\njulia> o = [bernoulliB(n) for n=0:5]; println(o)\nRational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]\n\njulia> o = bernoulliB(5; arr=true); println(o)\nRational{Int64}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]\n\njulia> o = bernoulliB(big(5); arr=true); println(o)\nRational{BigInt}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1]\n\njulia> bernoulliB(60)\nIOP capture: bernoulliB(60) converted to Rational{BigInt}\n-1215233140483755572040304994079820246041491//56786730\n\njulia> n = 60;\njulia> bernoulliB(n; msg=false) == bernoulliB(n; msg=false, arr=true)[end]             \ntrue\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Divisor-(common-denominator-of-Rationals)","page":"Library","title":"Divisor (common denominator of Rationals)","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"normalize_rationals(v::Vector{Rational{T}}) where {T<:Integer}\ndivisor(v::Vector{Rational{T}}) where {T<:Integer}\nnumerators(v::Vector{Rational{T}}) where {T<:Integer}","category":"page"},{"location":"man/library/#CamiMath.normalize_rationals-Union{Tuple{Array{Rational{T}, 1}}, Tuple{T}} where T<:Integer","page":"Library","title":"CamiMath.normalize_rationals","text":"normalize_rationals(v::Vector{Rational{T}}) where T<:Integer\n\nNumerators separated from divisor\n\nExample:\n\njulia> normalize_rationals([1//1, 1//2, 1//3])\n([6, 3, 2], 6)\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.divisor-Union{Tuple{Array{Rational{T}, 1}}, Tuple{T}} where T<:Integer","page":"Library","title":"CamiMath.divisor","text":"divisor(v::Vector{Rational{T}}) where {T<:Integer}\n\nGreatest common denominator of the set of rational numbers v\n\nExample:\n\njulia> divisor([1//1, 1//2, 1//3])\n6\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.numerators-Union{Tuple{Array{Rational{T}, 1}}, Tuple{T}} where T<:Integer","page":"Library","title":"CamiMath.numerators","text":"numerators(v::Vector{Rational{T}}) where {T<:Integer}\n\nNumerators for the standard devisor of the set of rational numbers v\n\nExample:\n\njulia> numerators([1//1, 1//2, 1//3])\n3-element Vector{Int64}:\n 6\n 3\n 2\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Factorial","page":"Library","title":"Factorial","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"bigfactorial(n::Integer; msg=true)","category":"page"},{"location":"man/library/#CamiMath.bigfactorial-Tuple{Integer}","page":"Library","title":"CamiMath.bigfactorial","text":"bigfactorial(n::Int [; msg=true])\n\nThe product of all positive integers less than or equal to n,\n\nn=n(n-1)(n-2)1\n\nIn addition 0=1 by definition. For negative integers the factorial is zero. \n\nmsg : integer-overflow protection (IOP) - warning on activation  (for n > 20) \n\nExamples:\n\njulia> bigfactorial(20) == factorial(20)\ntrue\n\njulia> bigfactorial(21)\nIOP capture: bigfactorial(21) converted to BigInt\n51090942171709440000\n\njulia> bigfactorial(21; msg=false)\n51090942171709440000\n\njulia> factorial(21)\nERROR: OverflowError: 21 is too large to look up in the table; consider using \n`factorial(big(21))` instead\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Faulhaber-polynomial","page":"Library","title":"Faulhaber polynomial","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"faulhaber_polynom(p::Integer; msg=true)\nfaulhaber_polynomial(n::Integer, p::Int; msg=true)\nfaulhaber_summation(n::Integer, p::Int)","category":"page"},{"location":"man/library/#CamiMath.faulhaber_polynom-Tuple{Integer}","page":"Library","title":"CamiMath.faulhaber_polynom","text":"faulhaber_polynom(p::Integer [; msg=true])\n\nVector representation of the coefficients of the faulhaber_polynomial  of degree p, \n\n   c=c_0 c_p\n\nwhere c_0=0  c_j=frac1pbinompp-jB_p-j, with j 1 p. The B_p-j are bernoulliB in the  even index convention (but with  B_1=+frac12 rather than -frac12).\n\nmsg : integer-overflow protection (IOP) - warning on activation \n\n(for p > 36)\n\nExample:\n\nfaulhaber_polynom(6)\n7-element Vector{Rational{Int64}}:\n  0//1\n  0//1\n -1//12\n  0//1\n  5//12\n  1//2\n  1//6\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.faulhaber_polynomial-Tuple{Integer, Int64}","page":"Library","title":"CamiMath.faulhaber_polynomial","text":"faulhaber_polynomial(n::Integer, p::Int [; msg=true])\n\nFaulhaber polynomial of degree p \n\n    F(np)=sum_j=0^pc_jn^j\n\nwhere n is a positive integer and the coefficients are contained in the  vector c=c_0 c_p given by faulhaber_polynom.\n\nmsg : integer-overflow protection (IOP) - warning on activation\n\nExamples:\n\njulia> faulhaber_polynomial(3, 6)\n276\n\njulia> faulhaber_polynomial(5, 30)\nIOP capture: faulhaber_polynomial(5, 30) autoconverted to Rational{BigInt}\n186552813930161650665\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.faulhaber_summation-Tuple{Integer, Int64}","page":"Library","title":"CamiMath.faulhaber_summation","text":"faulhaber_summation(n::Integer, p::Int [; msg=true])\n\nSum of the p^th power of the first n natural numbers\n\n    sum_k=1^nk^p=H_n-p=F(np+1)\n\nwhere H_n-p is a harmonicNumber  of power -p and F(np)  a faulhaber_polynomial of power p.\n\nmsg : integer-overflow protection (IOP) - warning on activation\n\nExamples:\n\njulia> faulhaber_summation(3,5)\n276\n\njulia> faulhaber_summation(3,60)\nIOP capture: faulhaber_polynom autoconverted to Rational{BigInt}\n42391158276369125018901280178\n\n\n\n\n\n","category":"method"},{"location":"man/library/#HarmonicNumber","page":"Library","title":"HarmonicNumber","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"harmonicNumber(n::Integer; arr=false, msg=true)\nharmonicNumber(n::Integer, p::Int; arr=false, msg=true)","category":"page"},{"location":"man/library/#CamiMath.harmonicNumber-Tuple{Integer}","page":"Library","title":"CamiMath.harmonicNumber","text":"harmonicNumber(n::Integer [, p=1 [; arr=false [], msg=true]]])\n\nSum of the p^th power of reciprocals of the first n positive integers,\n\n    H_np=sum_k=1^nfrac1k^p\n\narr : output in array format\nmsg : integer-overflow protection (IOP) - warning on activation \n\nExamples:\n\njulia> o = [harmonicNumber(n) for n=1:8]; println(o)\nRational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280]\n\njulia> @btime harmonicNumber(8; arr=true)\n(1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280)\n\njulia> @btime harmonicNumber(42)\n12309312989335019//2844937529085600\n\njulia> harmonicNumber(43)\nIOP capture: harmonicNumber(43, 1) converted to Rational{BigInt}\n532145396070491417//122332313750680800\n\njulia> harmonicNumber(12) == harmonicNumber(12, 1)\ntrue\n\nharmonicNumber(12, -3) == faulhaber_summation(12, 3)\n  true\n\njulia> o = [harmonicNumber(i, 5) for i=1:4]; println(o)\nRational{Int64}[1//1, 33//32, 8051//7776, 257875//248832]\n\njulia> o = harmonicNumber(4, 5; arr=true); println(o)\n(1//1, 33//32, 8051//7776, 257875//248832)\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.harmonicNumber-Tuple{Integer, Int64}","page":"Library","title":"CamiMath.harmonicNumber","text":"harmonicNumber(n::Integer [, p=1 [; arr=false [], msg=true]]])\n\nSum of the p^th power of reciprocals of the first n positive integers,\n\n    H_np=sum_k=1^nfrac1k^p\n\narr : output in array format\nmsg : integer-overflow protection (IOP) - warning on activation \n\nExamples:\n\njulia> o = [harmonicNumber(n) for n=1:8]; println(o)\nRational{Int64}[1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280]\n\njulia> @btime harmonicNumber(8; arr=true)\n(1//1, 3//2, 11//6, 25//12, 137//60, 49//20, 363//140, 761//280)\n\njulia> @btime harmonicNumber(42)\n12309312989335019//2844937529085600\n\njulia> harmonicNumber(43)\nIOP capture: harmonicNumber(43, 1) converted to Rational{BigInt}\n532145396070491417//122332313750680800\n\njulia> harmonicNumber(12) == harmonicNumber(12, 1)\ntrue\n\nharmonicNumber(12, -3) == faulhaber_summation(12, 3)\n  true\n\njulia> o = [harmonicNumber(i, 5) for i=1:4]; println(o)\nRational{Int64}[1//1, 33//32, 8051//7776, 257875//248832]\n\njulia> o = harmonicNumber(4, 5; arr=true); println(o)\n(1//1, 33//32, 8051//7776, 257875//248832)\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Fibonacci-number","page":"Library","title":"Fibonacci number","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"fibonacci(n::Integer; arr=false, msg=true)","category":"page"},{"location":"man/library/#CamiMath.fibonacci-Tuple{Integer}","page":"Library","title":"CamiMath.fibonacci","text":"fibonacci(n::Integer [[; arr=false], msg=true])\n\nThe sequence of integers,  F_0 F_nmax, in which each element is  the sum of the two preceding ones, \n\n    F_n = F_n-1+F_n-2\n\nwith F_1=1 and F_0=0. \n\narr : output full Pascal triangle\nmsg : integer-overflow protection (IOP) - warning on activation \n\n(for n > 92)\n\nExamples:\n\njulia> fibonacci(92)\n7540113804746346429\n\njulia> fibonacci(93)\nIOP capture: fibonaci(93) converted to BigInt\n12200160415121876738\n\njulia> o = fibonacci(10; arr=true); println(o)\n[1, 1, 2, 3, 5, 8, 13, 21, 34, 55]\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Integer-partitioning","page":"Library","title":"Integer partitioning","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"canonical_partitions(n::Int, m=0; header=true, reverse=true)\ninteger_partitions(n::Int, m=0; transpose=false, count=false)","category":"page"},{"location":"man/library/#CamiMath.canonical_partitions","page":"Library","title":"CamiMath.canonical_partitions","text":"canonical_partitions(n::Int, [[[m=0]; header=true,] reverse=true])\n\nCanonical partition of n in parts of maximum size m (m = 0 for any size)\n\nheader : unit partition included in output\n\nExamples:\n\njulia> canonical_partitions(6; header=true, reverse=false)\n6-element Vector{Vector{Int64}}:\n [6]\n [5, 1]\n [4, 2]\n [3, 3]\n [2, 2, 2]\n [1, 1, 1, 1, 1, 1]\n\njulia> canonical_partitions(6; header=true)\n6-element Vector{Vector{Int64}}:\n [1, 1, 1, 1, 1, 1]\n [2, 2, 2]\n [3, 3]\n [4, 2]\n [5, 1]\n [6]\n\njulia> canonical_partitions(6)\n6-element Vector{Vector{Int64}}:\n [1, 1, 1, 1, 1, 1]\n [2, 2, 2]\n [3, 3]\n [4, 2]\n [5, 1]\n [6]\n\njulia> o = canonical_partitions(9, 2); println(o)\n[2, 2, 2, 2, 1]\n\njulia> o = canonical_partitions(9, 3); println(o)\n[3, 3, 3]\n\n\n\n\n\n","category":"function"},{"location":"man/library/#CamiMath.integer_partitions","page":"Library","title":"CamiMath.integer_partitions","text":"integer_partitions(n [[[,m]; transpose=false], count=false])\n\ndefault                      : The integer partitions of n\n\ncount                        : The number of integer partitions of n\n\ntranspose = false (m > 0): partitions restricted to maximum part m             = true  (m > 0): partitions restricted to maximum length m`\n\ndefinitions:\n\nThe integer partition of the positive integer n is a nonincreasing sequence of positive integers p1, p2,... pk whose sum is n. The elements of the sequence are called the parts of the partition.\n\nExamples:\n\njulia> integer_partitions(7)\n15-element Vector{Vector{Int64}}:\n [1, 1, 1, 1, 1, 1, 1]\n [2, 2, 2, 1]\n [2, 2, 1, 1, 1]\n [2, 1, 1, 1, 1, 1]\n [3, 3, 1]\n [3, 2, 2]\n [3, 2, 1, 1]\n [3, 1, 1, 1, 1]\n [4, 3]\n [4, 2, 1]\n [4, 1, 1, 1]\n [5, 2]\n [5, 1, 1]\n [6, 1]\n [7]\n\njulia> integer_partitions(7; count=true)\n15\n\njulia> integer_partitions(7, 4; count=true)\n3\n\njulia> integer_partitions(7, 4)\n3-element Vector{Vector{Int64}}:\n [4, 3]\n [4, 2, 1]\n [4, 1, 1, 1]\n\njulia> integer_partitions(7, 4; transpose=true)\n3-element Vector{Vector{Int64}}:\n [2, 2, 2, 1]\n [3, 2, 1, 1]\n [4, 1, 1, 1]\n\n\n\n\n\n","category":"function"},{"location":"man/library/#Laguerre-polynomials","page":"Library","title":"Laguerre polynomials","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"laguerre_polynom(p::Integer; msg=true)\ngeneralized_laguerre_polynom(n::Int, α=0)","category":"page"},{"location":"man/library/#CamiMath.laguerre_polynom-Tuple{Integer}","page":"Library","title":"CamiMath.laguerre_polynom","text":"laguerre_polynom(n::Integer [; msg=true])\n\nThe coefficients of laguerreL for degree n, \n\n    v_n=c_0 c_1 cdots c_n\n\nwhere, with k=01n , \n\n    c_k = fracGamma(n+1)Gamma(k+1)frac(-1)^k(n-k)frac1k\n\nmsg : integer-overflow protection (IOP) - warning on activation \n\nExample:\n\njulia> laguerre_polynom(7)\n(1//1, -7//1, 21//2, -35//6, 35//24, -7//40, 7//720, -1//5040)\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.generalized_laguerre_polynom","page":"Library","title":"CamiMath.generalized_laguerre_polynom","text":"generalized_laguerre_polynom(n::Int [, α=0])\n\nThe coefficients of generalized_laguerreL for degree n and parameter α,\n\n    c(n α)k = fracGamma(α+n+1)Gamma(α+k+1)\n    frac(-1)^k(n-k)frac1k\n\nExample:\n\no = generalized_laguerre_polynom(8,3); println(o)\n    Rational{Int64}[165//1, -330//1, 231//1, -77//1, 55//4, -11//8, 11//144, -11//5040, 1//40320]\n\n\n\n\n\n","category":"function"},{"location":"man/library/#Pascal-triangle","page":"Library","title":"Pascal triangle","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"pascal_triangle(n::Integer; arr=false, msg=true)\npascal_next(a)","category":"page"},{"location":"man/library/#CamiMath.pascal_triangle-Tuple{Integer}","page":"Library","title":"CamiMath.pascal_triangle","text":"pascal_triangle(n::Integer [; arr=false [, msg=true]])\n\nRow n of Pascal triangle, r_n = binomn1cdots binomnn\n\narr : output full Pascal triangle\nmsg : integer-overflow protection (IOP) - warning on activation \n\n(for n > 10000)\n\nExamples:\n\njulia> [pascal_triangle(n) for n=0:5]\n6-element Vector{Vector{Int64}}:\n [1]\n [1, 1]\n [1, 2, 1]\n [1, 3, 3, 1]\n [1, 4, 6, 4, 1]\n [1, 5, 10, 10, 5, 1]\n\njulia> pascal_triangle(5; arr=true)\n5-element Vector{Vector{Int64}}:\n [1, 1]\n [1, 2, 1]\n [1, 3, 3, 1]\n [1, 4, 6, 4, 1]\n [1, 5, 10, 10, 5, 1]\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.pascal_next-Tuple{Any}","page":"Library","title":"CamiMath.pascal_next","text":"pascal_next(row)\n\nNext row of binomial coefficients of the Pascal triangle. \n\nExample:\n\njulia> pascal_next([1, 4, 6, 4, 1])\n6-element Vector{Int64}:\n  1\n  5\n 10\n 10\n  5\n  1\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Permutations","page":"Library","title":"Permutations","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"permutations_unique_count(p::Vector{Vector{Int}}, i::Int)","category":"page"},{"location":"man/library/#CamiMath.permutations_unique_count-Tuple{Vector{Vector{Int64}}, Int64}","page":"Library","title":"CamiMath.permutations_unique_count","text":"permutations_unique_count(p::Vector{Vector{Integer}}, i::Int)\n\nNumber of unique permutations of the subarray pi.\n\nExample:\n\np = [[1,2,3],[2,3,1,4,3]]\npermutations_unique_count(p,2)\n 60\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Pochhammer-product","page":"Library","title":"Pochhammer product","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"pochhammer(x::T, p::Int) where {T<:Real}","category":"page"},{"location":"man/library/#CamiMath.pochhammer-Union{Tuple{T}, Tuple{T, Int64}} where T<:Real","page":"Library","title":"CamiMath.pochhammer","text":"pochhammer(x::T, p::Int) where T<:Real\n\nPochhammer symbol (x)_p for non-negative integer p,\n\n(x)_p=begincases\n1  p=0\nx(x+1)(x+2)(x+p-1)  p0\nendcases\n\nNote that (x)_p=0 for x=0-1 -(p-1)\n\nExamples:\n\njulia> x = [-4,-3,-2,-1, 0, 1, 2 , 3, 4];\n\njulia> pochhammer.(x, 5) == [0, 0, 0, 0, 0, 120, 720, 2520, 6720]\ntrue\n\njulia> pochhammer.(x, 0) == [1, 1, 1, 1, 1, 1, 1, 1, 1]\ntrue\n\njulia> o = [pochhammer.([x for x=0:-1:-p],p) for p=0:5]\n6-element Vector{Vector{Int64}}:\n [1]\n [0, -1]\n [0, 0, 2]\n [0, 0, 0, -6]        \n [0, 0, 0, 0, 24]     \n [0, 0, 0, 0, 0, -120]\n\njulia>  o = [pochhammer.([x for x=0:p],p) for p=0:5]\n6-element Vector{Vector{Int64}}:\n [1]\n [0, 1]\n [0, 2, 6]\n [0, 6, 24, 60]\n [0, 24, 120, 360, 840]\n [0, 120, 720, 2520, 6720, 15120]\n\njulia> x = Rational{BigInt}(-1, 50);\n\njulia> pochhammer(x, 20)\n-21605762356630090481082546653745369902321614221999//9536743164062500000000000000000000\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Polynomials","page":"Library","title":"Polynomials","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"polynomial(coords, x::T; deriv=0) where {T<:Real}\npolynom_power(coords, power::Int)\npolynom_product(coords1, coords2)\npolynom_product_expansion(coords1, coords2, p::Int)","category":"page"},{"location":"man/library/#CamiMath.polynomial-Union{Tuple{T}, Tuple{Any, T}} where T<:Real","page":"Library","title":"CamiMath.polynomial","text":"polynomial(coords, x::T [; deriv=0]) where {T<:Real}\n\nPolynomial of degree d,\n\n    P(x)=c_0 + c_1 x +  + c_d x^d\n\nwhere coords = (c_0 c_d) are the coordinates defining the vector  representation of the polynomial in a vector space of dimension d+1.\n\nExamples:\n\njulia> coords = (1, 1, 1, 1, 1);\n           \njulia> P(x) = polynomial(coords,x);\n\njulia> P(1)\n5\n\njulia> polynomial(coords, 1; deriv=1)     # P′(1)\n10\n\njulia> polynomial(coords, 2; deriv=2)     # P″(1)\n20\n\njulia> polynomial(coords,x; deriv=-1)   # primitive (zero integration constant)\n137 // 60\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.polynom_power-Tuple{Any, Int64}","page":"Library","title":"CamiMath.polynom_power","text":"polynom_power(coords, p)\n\nCoordinates of the polynomial defined by coords raised to the power p, which define a polynomial in a vector space of dimension p d + 1, where d is the degree of the polynomial defined by coords.\n\nExamples:\n\njulia> coords = (1,1,1)    # coordinates of polynomial vector of degree d = 2\n(1, 1, 1)\n\njulia> coords = (1,1,1);\n\njulia> polynom_power(coords, 3)\n7-element Vector{Int64}:\n 1\n 3\n 6\n 7\n 6\n 3\n 1\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.polynom_product-Tuple{Any, Any}","page":"Library","title":"CamiMath.polynom_product","text":"polynom_product(coords1, coords2)\n\nCoordinate representation of the product of two polynomials, a = bcoords1 and  b = coords2 of degree m and n, which is a polynomial in a vector  space of dimension d=m+n+1,\n\n    P(x)=a_0b_0 + (a_0b_1 + b_0a_1)x +  + a_n b_m x^n+m\n\nExamples:\n\njulia> polynom_product((1, 1), (1, -1, 2))\n4-element Vector{Int64}:\n 1\n 0\n 1\n 2\n\njulia> polynom_product((1, 1), (1, -1.0, 2))\n4-element Vector{Real}:\n 1\n 0.0\n 1.0\n 2  \n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.polynom_product_expansion-Tuple{Any, Any, Int64}","page":"Library","title":"CamiMath.polynom_product_expansion","text":"polynom_product_expansion(coords1, coords2, p::Int)\n\nVector representation of the product of two polynomials, a = coords1  (of degree n) and b = coords2 (of degree m), with m  n truncated at the order p, which is a polynomial in a vector space of  dimension d=p+1\n\n\n\njulia> a = (1,-1,1);\n\njulia> b = (1,1,-1,1,1,1);\n\njulia> o = polynom_product(a, b); println(o)\n[1, 0, -1, 3, -1, 1, 0, 1]\n \njulia> o = polynom_product_expansion(a, b, 4); println(o)\n[1, 0, -1, 3, -1] \n\n\n\n\n\n","category":"method"},{"location":"man/library/#Triangle-relations","page":"Library","title":"Triangle relations","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"istriangle(a::Real, b::Real, c::Real)\ntriangle_coefficient(a::Real, b::Real, c::Real)","category":"page"},{"location":"man/library/#CamiMath.istriangle-Tuple{Real, Real, Real}","page":"Library","title":"CamiMath.istriangle","text":"istriangle(a::Real, b::Real, c::Real)\n\nTriangle condition for a triangle of sides a, b and c.\n\nExample:\n\njulia> istriangle(3, 4, 5)\ntrue\n\njulia> istriangle(1//2, 1, 1.5)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.triangle_coefficient-Tuple{Real, Real, Real}","page":"Library","title":"CamiMath.triangle_coefficient","text":"triangle_coefficient(a::Real, b::Real, c::Real)\n\nTriangle coefficient for a triangle of sides a, b and c,\n\n    Delta(abc)equivfrac(a+b-c)(b+c-a)(c+a-b)(a+b+c+1)\n\nTriangle condition satisfied for Delta  0\n\nExamples:\n\njulia> triangle_coefficient(3, 4, 5)\n1//180180\n\njulia> triangle_coefficient(1//2, 1, 1.5)\n1//12\n\n\n\n\n\n","category":"method"},{"location":"man/library/#Truncated-exponentials","page":"Library","title":"Truncated exponentials","text":"","category":"section"},{"location":"man/library/","page":"Library","title":"Library","text":"texp(x::T, a::T, p::Int) where {T<:Real}\nlog10_characteristic_power(x)\nlog10_mantissa(x)","category":"page"},{"location":"man/library/#CamiMath.texp-Union{Tuple{T}, Tuple{T, T, Int64}} where T<:Real","page":"Library","title":"CamiMath.texp","text":"texp(x::T, a::T, p::Int) where T <: Real\n\nTruncated exponential: Taylor expansion of exp(x) about x = a  up to order p`,\n\n    mathsftexp(xap) = 1+(x-a)+frac12(x-a)^2++frac1p(x-a)^p\n\nExamples:\n\njulia> texp(1.0, 0.0, 5)\n2.7166666666666663\n\njulia> texp(1, 0, 5)\n163//60\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.log10_characteristic_power-Tuple{Any}","page":"Library","title":"CamiMath.log10_characteristic_power","text":"log10_characteristic_power(x)\n\ncharacteristic power-of-10 of the number x\n\nExamples:\n\nlog10_characteristic_power.([3,30,300])\n3-element Vector{Int64}:\n 0\n 1\n 2\n\n\n\n\n\n","category":"method"},{"location":"man/library/#CamiMath.log10_mantissa-Tuple{Any}","page":"Library","title":"CamiMath.log10_mantissa","text":"log10_mantissa(x)\n\nlog10 mantissa of the number x\n\nExamples:\n\nlog10_mantissa.([3,30,300])\n3-element Vector{Float64}:\n 0.47712125471966244\n 0.4771212547196624\n 0.4771212547196626\n\n\n\n\n\n","category":"method"},{"location":"#CamiMath.jl","page":"Home","title":"CamiMath.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Mathematics library with integer-overload protection (IOP)","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Library","page":"Home","title":"Library","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"man/library.md\"]","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"man/index.md\"]","category":"page"},{"location":"man/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"man/","page":"Index","title":"Index","text":"","category":"page"}]
}
