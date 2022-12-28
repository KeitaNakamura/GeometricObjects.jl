@testset "Geometry" begin
    @testset "misc" begin
        x = Line(Vec(0.0,0.0) => Vec(1.0,0.0))
        @test coordinates(x, 1) == Vec(0,0)
        @test coordinates(x, 2) == Vec(1,0)
    end
end
