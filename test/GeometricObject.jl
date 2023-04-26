@testset "GeometricObject" begin
    @testset "Geometry" begin
        line = Line(Vec(0.0,0.0) => Vec(1.0,0.0))
        obj = @inferred GeometricObject(line)
        # geometry
        @test (@inferred geometry(obj)) === line
        # coordinates
        @test (@inferred coordinates(obj)) === coordinates(line)
        @test (@inferred coordinates(obj, 1)) === coordinates(line, 1)
        @test (@inferred coordinates(obj, 2)) === coordinates(line, 2)
        @test_throws Exception coordinates(obj, 3) # out of range
        # quaternion
        @test (@inferred quaternion(obj)) === quaternion(line)
        # attitude
        @test (@inferred attitude(obj)) === attitude(line)
    end
    @testset "apply_force!" begin
        ## 2D
        line = GeometricObject(Line(Vec(0.3,0.4) => Vec(0.5,0.4)))
        v = coordinates(line, 2) - coordinates(line, 1)
        l² = v ⋅ v
        I = l² / 12 * line.m
        update_geometry!(apply_force!(line, Vec(-0.8,-0.8), I, 1.0), 1.0)
        @test moment_of_inertia(geometry(line)) ≈ I
        @test norm(coordinates(line, 2) - coordinates(line, 1)) ≈ norm(v)
        @test (x = coordinates(line, 2) - coordinates(line, 1); atan(x[2]/x[1])) ≈ 1
        ## 3D
        line = GeometricObject(Line(Vec(0.0, 0.0, 0.0) => Vec(1.0, 0.0, 0.0)))
        v = coordinates(line, 2) - coordinates(line, 1)
        l² = v ⋅ v
        I = l² / 12 * line.m
        update_geometry!(apply_force!(line, Vec(0.5,0.5,0.5), Vec(0.0,0.0,I), 1.0), 1.0)
        @test centroid(geometry(line)) == [1.0, 0.5, 0.5]
        @test norm(coordinates(line, 2) - coordinates(line, 1)) ≈ norm(v)
        @test (x = coordinates(line, 2) - coordinates(line, 1); atan(x[2]/x[1])) ≈ 1
    end
end
