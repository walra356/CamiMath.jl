using CamMath
using Test

@testset "CamMath.jl" begin
    # Write your tests here.
    @test fA(1) == 1
    @test fA(2) == 4
    @test fA(3) == 9
    @test fA(4) == 16
    @test fA(5) == 25
end
