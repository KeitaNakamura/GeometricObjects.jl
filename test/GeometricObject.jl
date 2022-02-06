@testset "GeometricObject" begin
    @testset "misc" begin
        x = GeometricObject(Line(Vec(0.0,0.0) => Vec(1.0,0.0)))
        @test eachindex(x) == Base.OneTo(2)
        @test firstindex(x) == 1
        @test lastindex(x) == 2
        @test x[begin] == Vec(0,0)
        @test x[end] == Vec(1,0)
        # @test (@inferred collect(x))::Vector{Vec{2, Float64}} ≈ [Vec(0,0), Vec(1,0)]
    end
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
        v = line[2] - line[1]
        l² = v ⋅ v
        I = l² / 12 * line.m
        for (F, τ) in ((Vec(-0.8,-0.8), Vec(0.0,0.0,I)), ([-0.8,-0.8], [0.0,0.0,I]))
            line = GeometricObject(Line(Vec(0.3,0.4) => Vec(0.5,0.4)))
            GeometricObjects.update!(line, F, τ, 1.0)
            @test moment_of_inertia(line)[3,3] ≈ I
            @test norm(line[2] - line[1]) ≈ norm(v)
            @test (x = line[2] - line[1]; atan(x[2]/x[1])) ≈ 1
        end
        ## 3D
        line = GeometricObject(Line(Vec(0.0, 0.0, 0.0) => Vec(1.0, 0.0, 0.0)))
        v = line[2] - line[1]
        l² = v ⋅ v
        I = l² / 12 * line.m
        for (F, τ) in ((Vec(0.5,0.5,0.5), Vec(0.0,0.0,I)), ([0.5,0.5,0.5], [0.0,0.0,I]))
            line = GeometricObject(Line(Vec(0.0, 0.0, 0.0) => Vec(1.0, 0.0, 0.0)))
            GeometricObjects.update!(line, Vec(0.5,0.5,0.5), Vec(0.0,0.0,I), 1.0)
            @test centroid(line) == [1.0, 0.5, 0.5]
            @test norm(line[2] - line[1]) ≈ norm(v)
            @test (x = line[2] - line[1]; atan(x[2]/x[1])) ≈ 1
        end
    end
end
