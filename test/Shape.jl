@testset "Shape" begin
    @testset "misc" begin
        x = Line(Vec(0.0,0.0) => Vec(1.0,0.0))
        @test eachindex(x) == Base.OneTo(2)
        @test firstindex(x) == 1
        @test lastindex(x) == 2
        @test x[begin] == Vec(0,0)
        @test x[end] == Vec(1,0)
        # @test (@inferred collect(x))::Vector{Vec{2, Float64}} ≈ [Vec(0,0), Vec(1,0)]
    end
end
