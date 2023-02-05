using CamMath
using Test

@testset "CamMath.jl" begin

    @test bernoulliB_array(10) == [1 // 1, -1 // 2, 1 // 6, 0 // 1, -1 // 30, 0 // 1, 1 // 42, 0 // 1, -1 // 30, 0 // 1, 5 // 66]
    @test bernoulliB_array(big(8)) == Rational{BigInt}[1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30]
    @test bernoulliB(60; msg=false) == -1215233140483755572040304994079820246041491 // 56786730
    @test bernoulliB(60; msg=false) == bernoulliB_array(60; msg=false)[end]
    @test bigfactorial(21; msg=false) == 51090942171709440000
    @test sum([sum(faulhaber_polynom(p; msg=false)) for p = 1:40]) == 40 // 1
    @test sum([sum(faulhaber_polynom(big(p); msg=false)) for p = 1:40]) == 40 // 1
    #   ............................................................................
    @test faulhaber_polynomial(3, 6) == 276
    @test faulhaber_polynomial(5, 30; msg=false) == 186552813930161650665
    @test faulhaber_polynomial(5, 36; msg=false) == (2911563687325667231369625)
    @test faulhaber_summation(3, 5) == 276
    #   ............................................................................
    @test harmonicNumber(3, -5) == 276
    @test harmonicNumber(46) == (5943339269060627227 // 1345655451257488800)
    @test harmonicNumber(47; msg=false) == (280682601097106968469 // 63245806209101973600)
    @test harmonicNumber(24, 2) == (187700554334941861 // 117011293467045120)
    @test harmonicNumber(25, 2; msg=false) == (23485971550561141649 // 14626411683380640000)
    @test harmonicNumber(2, 10) == (1025 // 1024)
    @test harmonicNumber(2, 11; msg=false) == (2049 // 2048)
    @test typeof(harmonicNumber(8)) == Rational{Int}
    @test typeof(harmonicNumber(big(8))) == Rational{BigInt}
    @test typeof(harmonicNumber(12, 3)) == Rational{Int}
    @test typeof(harmonicNumber(big(12), 3)) == Rational{BigInt}
    #   ............................................................................
    @test harmonicNumber_array(8) == [1 // 1, 3 // 2, 11 // 6, 25 // 12, 137 // 60, 49 // 20, 363 // 140, 761 // 280]
    @test harmonicNumber_array(big(8)) == [1 // 1, 3 // 2, 11 // 6, 25 // 12, 137 // 60, 49 // 20, 363 // 140, 761 // 280]
    @test pascal_triangle(5) == [[1], [1, 1], [1, 2, 1], [1, 3, 3, 1], [1, 4, 6, 4, 1], [1, 5, 10, 10, 5, 1]]
    @test pascal_next([1, 4, 6, 4, 1]) == [1, 5, 10, 10, 5, 1]

end
