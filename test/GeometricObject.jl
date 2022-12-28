@testset "GeometricObject" begin
    @testset "compute_force_moment" begin
        circle = GeometricObject(Circle(Vec(0.0,0.0), 2.0))
        F, M = GeometricObjects.compute_force_moment(circle, [Vec(0.0,1.0), Vec(0.0,1.0)], [Vec(1.0,0.0), Vec(-0.5,0.0)])
        @test F ≈ Vec(0.0,2.0)
        @test M ≈ Vec(0.0,0.0,0.5)
        F, M = GeometricObjects.compute_force_moment(circle, Vec{2,Float64}[], Vec{2,Float64}[])
        @test F ≈ Vec(0.0,0.0)
        @test M ≈ Vec(0.0,0.0,0.0)
    end
    # check common methods for GeometricObject
    @testset "update!" begin
        ## 2D
        line = GeometricObject(Line(Vec(0.3,0.4) => Vec(0.5,0.4)))
        v = coordinates(line, 2) - coordinates(line, 1)
        l² = v ⋅ v
        I = l² / 12 * line.m
        GeometricObjects.update!(line, Vec(-0.8,-0.8), Vec(0.0,0.0,I), 1.0)
        @test moment_of_inertia(geometry(line))[3,3] ≈ I
        @test norm(coordinates(line, 2) - coordinates(line, 1)) ≈ norm(v)
        @test (x = coordinates(line, 2) - coordinates(line, 1); atan(x[2]/x[1])) ≈ 1
        ## 3D
        line = GeometricObject(Line(Vec(0.0, 0.0, 0.0) => Vec(1.0, 0.0, 0.0)))
        v = coordinates(line, 2) - coordinates(line, 1)
        l² = v ⋅ v
        I = l² / 12 * line.m
        GeometricObjects.update!(line, Vec(0.5,0.5,0.5), Vec(0.0,0.0,I), 1.0)
        @test centroid(geometry(line)) == [1.0, 0.5, 0.5]
        @test norm(coordinates(line, 2) - coordinates(line, 1)) ≈ norm(v)
        @test (x = coordinates(line, 2) - coordinates(line, 1); atan(x[2]/x[1])) ≈ 1
    end
end
