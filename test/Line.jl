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
    @test centroid(line) == (line[1] + line[2]) / 2
    @test GeometricObjects.normalunit(line) ≈ [0.0,-1.0]
    # 1
    @test distance(line, Vec(0.2,0.6)) ≈ [0.0,-0.2]
    @test distance(line, Vec(0.2,0.6), r) ≈ a - [0.2,0.6]
    @test GeometricObjects.distance_from_outside(line, Vec(0.2,0.6), r) === nothing
    @test GeometricObjects.perpendicularfoot(line, Vec(0.2,0.6)) ≈ [0.2,0.4]
    # 2
    @test distance(line, Vec(0.4,0.6)) ≈ [0.0,-0.2]
    @test distance(line, Vec(0.4,0.6), r) ≈ [0.0,-0.2]
    @test GeometricObjects.distance_from_outside(line, Vec(0.4,0.6), r) === nothing
    @test GeometricObjects.perpendicularfoot(line, Vec(0.4,0.6)) ≈ [0.4,0.4]
    # 3
    @test distance(line, Vec(0.6,0.6)) ≈ [0.0,-0.2]
    @test distance(line, Vec(0.6,0.6), r) ≈ b - [0.6,0.6]
    @test GeometricObjects.distance_from_outside(line, Vec(0.6,0.6), r) === nothing
    @test GeometricObjects.perpendicularfoot(line, Vec(0.6,0.6)) ≈ [0.6,0.4]
    # 4
    @test distance(line, Vec(0.2,0.2)) ≈ [0.0,0.2]
    @test distance(line, Vec(0.2,0.2), r) ≈ a - [0.2,0.2]
    @test GeometricObjects.distance_from_outside(line, Vec(0.2,0.2), r) ≈ a - [0.2,0.2]
    @test GeometricObjects.perpendicularfoot(line, Vec(0.2,0.2)) ≈ [0.2,0.4]
    # 5
    @test distance(line, Vec(0.4,0.2)) ≈ [0.0,0.2]
    @test distance(line, Vec(0.4,0.2), r) ≈ [0.0,0.2]
    @test GeometricObjects.distance_from_outside(line, Vec(0.4,0.2), r) ≈ [0.4,0.4] - [0.4,0.2]
    @test GeometricObjects.perpendicularfoot(line, Vec(0.4,0.2)) ≈ [0.4,0.4]
    # 6
    @test distance(line, Vec(0.6,0.2)) ≈ [0.0,0.2]
    @test distance(line, Vec(0.6,0.2), r) ≈ b - [0.6,0.2]
    @test GeometricObjects.distance_from_outside(line, Vec(0.6,0.2), r) ≈ [0.5,0.4] - [0.6,0.2]
    @test GeometricObjects.perpendicularfoot(line, Vec(0.6,0.2)) ≈ [0.6,0.4]

    # enlarge
    line = Line(Vec(0.0, 0.0), (Vec(1.0, 0.0)))
    @test (@inferred enlarge(line, 1.1))::Line ≈ [[-0.05, 0.0], [1.05, 0.0]]
end
