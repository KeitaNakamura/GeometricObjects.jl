@testset "Shape" begin
    @testset "misc" begin
        x = Line(Vec(0.0,0.0) => Vec(1.0,0.0))
        @test eachindex(x) == Base.OneTo(2)
        # @test (@inferred collect(x))::Vector{Vec{2, Float64}} ≈ [Vec(0,0), Vec(1,0)]
    end
end
