@testset "GeometricObject" begin
    # check common methods for GeometricObject
    @testset "Line" begin
        ## 2D
        line = GeometricObject(Line(Vec(0.3,0.4) => Vec(0.5,0.4)))
        v = line[2] - line[1]
        l² = v ⋅ v
        I = l² / 12 * line.m
        GeometricObjects.update!(line, Vec(-0.8,-0.8), Vec(0.0,0.0,I), 1.0)
        @test moment_of_inertia(line)[3,3] ≈ I
        @test norm(line[2] - line[1]) ≈ norm(v)
        @test (x = line[2] - line[1]; atan(x[2]/x[1])) ≈ 1
        ## 3D
        line = GeometricObject(Line(Vec(0.0, 0.0, 0.0) => Vec(1.0, 0.0, 0.0)))
        v = line[2] - line[1]
        l² = v ⋅ v
        I = l² / 12 * line.m
        GeometricObjects.update!(line, Vec(0.5,0.5,0.5), Vec(0.0,0.0,I), 1.0)
        @test centroid(line) == [1.0, 0.5, 0.5]
        @test norm(line[2] - line[1]) ≈ norm(v)
        @test (x = line[2] - line[1]; atan(x[2]/x[1])) ≈ 1
    end
end
