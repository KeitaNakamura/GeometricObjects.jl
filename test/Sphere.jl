@testset "Circle" begin
    circle = Circle(Vec(1.0,1.0), 2.0)
    @test area(circle) ≈ π*2^2
end
