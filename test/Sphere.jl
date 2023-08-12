@testset "Sphere" begin
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

    # rotate! with origin
    sphere = Sphere(Vec(1.0,1.0,1.0), 2.0)
    rotate!(sphere, Vec(0,0,π), Vec(0.0,0.0,0.0))
    @test centroid(sphere) ≈ [-1.0,-1.0,1.0]
end

@testset "Circle" begin
    circle = Circle(Vec(1.0,1.0), 2.0)
    @test area(circle) ≈ π*2^2

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
