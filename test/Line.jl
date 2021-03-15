@testset "Line" begin
    #
    #     ^
    #     |
    # 0.6 +   1       2       3
    #     |
    # 0.4 +       a-------b
    #     |
    # 0.2 +   4       5       6
    #     |
    # ----+---+-------+-------+------>
    #     |  0.2     0.4     0.6
    #     |
    a = Vec(0.3,0.4)
    b = Vec(0.5,0.4)
    r = 0.3
    line = @inferred Line(a, b)
    @test (@inferred Line(a => b)) == line
    line.m = 2.0
    @test centroid(line) == (line[1] + line[2]) / 2
    @test normalunit(line) ≈ [0.0,-1.0]
    # 1
    @test distance(line, Vec(0.2,0.6)) ≈ [0.0,0.4] - [0.0,0.6]
    @test distance(line, Vec(0.2,0.6), r) === nothing
    @test perpendicularfoot(line, Vec(0.2,0.6)) ≈ [0.2,0.4]
    # 2
    @test distance(line, Vec(0.4,0.6)) ≈ [0.0,0.4] - [0.0,0.6]
    @test distance(line, Vec(0.4,0.6), r) === nothing
    @test perpendicularfoot(line, Vec(0.4,0.6)) ≈ [0.4,0.4]
    # 3
    @test distance(line, Vec(0.6,0.6)) ≈ [0.0,0.4] - [0.0,0.6]
    @test distance(line, Vec(0.6,0.6), r) === nothing
    @test perpendicularfoot(line, Vec(0.6,0.6)) ≈ [0.6,0.4]
    # 4
    @test distance(line, Vec(0.2,0.2)) ≈ [0.0,0.4] - [0.0,0.2]
    @test distance(line, Vec(0.2,0.2), r) ≈ [0.3,0.4] - [0.2,0.2]
    @test perpendicularfoot(line, Vec(0.2,0.2)) ≈ [0.2,0.4]
    # 5
    @test distance(line, Vec(0.4,0.2)) ≈ [0.0,0.4] - [0.0,0.2]
    @test distance(line, Vec(0.4,0.2), r) ≈ [0.4,0.4] - [0.4,0.2]
    @test perpendicularfoot(line, Vec(0.4,0.2)) ≈ [0.4,0.4]
    # 6
    @test distance(line, Vec(0.6,0.2)) ≈ [0.0,0.4] - [0.0,0.2]
    @test distance(line, Vec(0.6,0.2), r) ≈ [0.5,0.4] - [0.6,0.2]
    @test perpendicularfoot(line, Vec(0.6,0.2)) ≈ [0.6,0.4]

    # check update! function
    v = b - a
    l² = v ⋅ v
    I = l² / 12 * line.m
    GeometricObjects.update!(line, Vec(-0.8,-0.8), Vec(0.0,0.0,I), 1.0)
    @test moment_of_inertia(line)[3,3] ≈ I
    @test norm(line[2] - line[1]) ≈ norm(b - a)
    @test (x = line[2] - line[1]; atan(x[2]/x[1])) ≈ 1
end
