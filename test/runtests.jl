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
    @test faulhaber_polynomial(3, 6) == 276
    @test faulhaber_polynomial(5, 30; msg=false) == 186552813930161650665
    @test faulhaber_summation(3, 5) == 276
    @test harmonicNumber(3, -5) == 276
    @test harmonicNumber(1) == 1 // 1
    @test harmonicNumber(60) == 15117092380124150817026911 // 3230237388259077233637600
    @test harmonicNumber(12, 3) == 25535765062457 // 21300003648000
    @test harmonicNumber_array(9) == [1 // 1, 3 // 2, 11 // 6, 25 // 12, 137 // 60, 49 // 20, 363 // 140, 761 // 280, 7129 // 2520]
    @test pascal_triangle(5) == [[1], [1, 1], [1, 2, 1], [1, 3, 3, 1], [1, 4, 6, 4, 1], [1, 5, 10, 10, 5, 1]]
    @test pascal_next([1, 4, 6, 4, 1]) == [1, 5, 10, 10, 5, 1]

end
