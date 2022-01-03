@testset "Polyline" begin
    poly = Polyline(Vec(0.0, 0.0), Vec(1.0, 10.0), Vec(0.0, 20.0))
    @test distance(poly, Vec(0.0, -1.0), 2.0, include_tips = true)  == poly[1] - Vec(0.0, -1.0)
    @test distance(poly, Vec(0.0, -1.0), 2.0, include_tips = false) == nothing
    @test distance(poly, Vec(1.0, 0.0), 2.0, include_tips = true)  == distance(GeometricObjects.getline(poly, 1), Vec(1.0, 0.0), 2.0)
    @test distance(poly, Vec(1.0, 0.0), 2.0, include_tips = false) == distance(GeometricObjects.getline(poly, 1), Vec(1.0, 0.0), 2.0)
    @test distance(poly, Vec(2.0, 10.0), 2.0, include_tips = true)  == Vec(-1.0, 0.0)
    @test distance(poly, Vec(2.0, 10.0), 2.0, include_tips = false) == Vec(-1.0, 0.0)
    @test distance(poly, Vec(-1.0, 20.0), 2.0, include_tips = true)  == poly[3] - Vec(-1.0, 20.0)
    @test distance(poly, Vec(-1.0, 20.0), 2.0, include_tips = false) == nothing
    @test distance(poly, Vec(0.0, 21.0), 2.0, include_tips = true)  == poly[3] - Vec(0.0, 21.0)
    @test distance(poly, Vec(0.0, 21.0), 2.0, include_tips = false) == nothing

    # intersect
    poly = Polyline(Vec(0.0,1.0), Vec(1.0,0.0), Vec(1.0,1.0))
    line = Line(Vec(0.0,0.0), Vec(0.2,0.2));
    @test intersect(poly, line, extended = true) ≈ [0.5,0.5]
    @test intersect(poly, line, extended = false) === nothing
    line = Line(Vec(2.0,-0.2), Vec(1.8,-0.2));
    @test intersect(poly, line, extended = true) === nothing
    @test intersect(poly, line, extended = false) === nothing
end
