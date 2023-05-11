@testset "Sphere" begin
    sphere = Sphere(Vec(1.0,1.0,1.0), 2.0)

    # enlarge
    R = 1.1
    sphere′ = (@inferred enlarge(sphere, R))::Sphere
    @test coordinates(sphere) == coordinates(sphere′)
    @test centroid(sphere) == centroid(sphere′)
    @test R * radius(sphere) == radius(sphere′)

    # centered
    @test (@inferred GeometricObjects.centered(sphere))::Sphere |> coordinates ≈ [[0.0, 0.0, 0.0]]

    # distance
    sphere = Sphere(Vec(3.0,3.0,3.0), 2√3)
    @test distance(sphere, Vec(0,0,0), √3+1) ≈ [1,1,1]
    @test distance(sphere, Vec(0,0,0), √3-1) === nothing
    @test distance(sphere, Vec(2,2,2), √3+1) ≈ [-1,-1,-1]
    @test distance(sphere, Vec(2,2,2), √3-1) ≈ [-1,-1,-1]
    @test distance(sphere, Vec(0,0,0), √3+1; inverse=true) ≈ [1,1,1]
    @test distance(sphere, Vec(0,0,0), √3-1; inverse=true) ≈ [1,1,1]
    @test distance(sphere, Vec(2,2,2), √3+1; inverse=true) ≈ [-1,-1,-1]
    @test distance(sphere, Vec(2,2,2), √3-1; inverse=true) === nothing
end

@testset "Circle" begin
    circle = Circle(Vec(1.0,1.0), 2.0)
    @test area(circle) ≈ π*2^2

    # enlarge
    R = 1.1
    circle′ = (@inferred enlarge(circle, R))::Circle
    @test coordinates(circle) == coordinates(circle′)
    @test centroid(circle) == centroid(circle′)
    @test R * radius(circle) == radius(circle′)

    # centered
    @test (@inferred GeometricObjects.centered(circle))::Circle |> coordinates ≈ [[0.0, 0.0]]

    # distance
    sphere = Circle(Vec(3.0,3.0), 2√2)
    @test distance(sphere, Vec(0,0), √2+1) ≈ [1,1]
    @test distance(sphere, Vec(0,0), √2-1) === nothing
    @test distance(sphere, Vec(2,2), √2+1) ≈ [-1,-1]
    @test distance(sphere, Vec(2,2), √2-1) ≈ [-1,-1]
    @test distance(sphere, Vec(0,0), √2+1; inverse=true) ≈ [1,1]
    @test distance(sphere, Vec(0,0), √2-1; inverse=true) ≈ [1,1]
    @test distance(sphere, Vec(2,2), √2+1; inverse=true) ≈ [-1,-1]
    @test distance(sphere, Vec(2,2), √2-1; inverse=true) === nothing
end
