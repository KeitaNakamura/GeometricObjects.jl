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
    c = coordinates(line)
    @test coordinates(@inferred Line(a => b)) == coordinates(line)
    @test centroid(line) == (c[1] + c[2]) / 2
    @test GeometricObjects.normalunit(line) ≈ [0.0,-1.0]
    # 1
    @test distance(line, Vec(0.2,0.6)) ≈ [0.0,-0.2]
    @test distance(line, Vec(0.2,0.6), r) === nothing
    @test GeometricObjects.perpendicularfoot(line, Vec(0.2,0.6)) ≈ [0.2,0.4]
    # 2
    @test distance(line, Vec(0.4,0.6)) ≈ [0.0,-0.2]
    @test distance(line, Vec(0.4,0.6), r) ≈ [0.0,-0.2]
    @test GeometricObjects.perpendicularfoot(line, Vec(0.4,0.6)) ≈ [0.4,0.4]
    # 3
    @test distance(line, Vec(0.6,0.6)) ≈ [0.0,-0.2]
    @test distance(line, Vec(0.6,0.6), r) === nothing
    @test GeometricObjects.perpendicularfoot(line, Vec(0.6,0.6)) ≈ [0.6,0.4]
    # 4
    @test distance(line, Vec(0.2,0.2)) ≈ [0.0,0.2]
    @test distance(line, Vec(0.2,0.2), r) === nothing
    @test GeometricObjects.perpendicularfoot(line, Vec(0.2,0.2)) ≈ [0.2,0.4]
    # 5
    @test distance(line, Vec(0.4,0.2)) ≈ [0.0,0.2]
    @test distance(line, Vec(0.4,0.2), r) ≈ [0.0,0.2]
    @test GeometricObjects.perpendicularfoot(line, Vec(0.4,0.2)) ≈ [0.4,0.4]
    # 6
    @test distance(line, Vec(0.6,0.2)) ≈ [0.0,0.2]
    @test distance(line, Vec(0.6,0.2), r) === nothing
    @test GeometricObjects.perpendicularfoot(line, Vec(0.6,0.2)) ≈ [0.6,0.4]

    # norm
    line = Line(Vec(0.0,0.0), Vec(1.0,1.0))
    c = coordinates(line)
    @test (@inferred norm(line)) ≈ norm(c[1] - c[2])

    # intersect
    ## 2D
    line1 = Line(Vec(0.0,0.0), Vec(0.2,0.2));
    line2 = Line(Vec(0.0,1.0), Vec(1.0,0.0));
    @test intersect(line1, line2, extended = (true, true))  ≈ [0.5,0.5]
    @test intersect(line1, line2, extended = (true, false)) ≈ [0.5,0.5]
    @test intersect(line1, line2, extended = (false, true))  === nothing
    @test intersect(line1, line2, extended = (false, false)) === nothing
    line1 = Line(Vec(0.0,0.0), Vec(0.2,0.2));
    line2 = Line(Vec(0.0,1.0), Vec(0.2,0.8));
    @test intersect(line1, line2, extended = (true, true))  ≈ [0.5,0.5]
    @test intersect(line1, line2, extended = (true, false)) === nothing
    @test intersect(line1, line2, extended = (false, true))  === nothing
    @test intersect(line1, line2, extended = (false, false)) === nothing
    ## 3D
    line1 = Line(Vec(0.0,0.0,0.0), Vec(0.2,0.2,0.2));
    line2 = Line(Vec(0.0,1.0,0.0), Vec(1.0,0.0,1.0));
    @test intersect(line1, line2, extended = (true, true))  ≈ [0.5,0.5,0.5]
    @test intersect(line1, line2, extended = (true, false)) ≈ [0.5,0.5,0.5]
    @test intersect(line1, line2, extended = (false, true))  === nothing
    @test intersect(line1, line2, extended = (false, false)) === nothing
end
