@testset "Sphere" begin
    sphere = Sphere(Vec(1.0,1.0,1.0), 2.0)

    # enlarge
    R = 1.1
    sphere′ = (@inferred enlarge(sphere, R))::Sphere
    @test sphere == sphere′
    @test centroid(sphere) == centroid(sphere′)
    @test R * radius(sphere) == radius(sphere′)
end

@testset "Circle" begin
    circle = Circle(Vec(1.0,1.0), 2.0)
    @test area(circle) ≈ π*2^2

    # enlarge
    R = 1.1
    circle′ = (@inferred enlarge(circle, R))::Circle
    @test circle == circle′
    @test centroid(circle) == centroid(circle′)
    @test R * radius(circle) == radius(circle′)
end
