using CamMath
using Test

@testset "CamMath.jl" begin

    @test bernoulliB_array(10) == [1 // 1, -1 // 2, 1 // 6, 0 // 1, -1 // 30, 0 // 1, 1 // 42, 0 // 1, -1 // 30, 0 // 1, 5 // 66] 
    @test bernoulliB(0) == 1 // 1
    @test bernoulliB(1) == -1 // 2
    @test bernoulliB(60; msg=false) == -1215233140483755572040304994079820246041491 // 56786730
    @test (bernoulliB(60; msg=false) == bernoulliB_array(60; msg=false)[end]) == true

end
