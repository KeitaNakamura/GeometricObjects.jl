@testset "Polygon" begin
    # `in`
    poly = Polygon([Vec(0.1,0.1), Vec(0.3,0.1), Vec(0.4,0.0), Vec(0.5,0.1), Vec(0.7,0.1),
                    Vec(0.6,0.4), Vec(0.5,0.4), Vec(0.4,0.3), Vec(0.3,0.3), Vec(0.2,0.4)])
    # Vec(0.0,0.0)
    @test in(Vec(0.0,0.0), poly; include_bounds = true)  == false
    @test in(Vec(0.0,0.0), poly; include_bounds = false) == false
    # Vec(0.0,0.1)
    @test in(Vec(0.0,0.1), poly; include_bounds = true)  == false
    @test in(Vec(0.0,0.1), poly; include_bounds = false) == false
    # Vec(0.2,0.1)
    @test in(Vec(0.2,0.1), poly; include_bounds = true)  == true
    @test in(Vec(0.2,0.1), poly; include_bounds = false) == false
    # Vec(0.4,0.1)
    @test in(Vec(0.4,0.1), poly; include_bounds = true)  == true
    @test in(Vec(0.4,0.1), poly; include_bounds = false) == true
    # Vec(0.6,0.1)
    @test in(Vec(0.6,0.1), poly; include_bounds = true)  == true
    @test in(Vec(0.6,0.1), poly; include_bounds = false) == false
    # Vec(0.1,0.2)
    @test in(Vec(0.1,0.2), poly; include_bounds = true)  == false
    @test in(Vec(0.1,0.2), poly; include_bounds = false) == false
    # Vec(0.3,0.2)
    @test in(Vec(0.3,0.2), poly; include_bounds = true)  == true
    @test in(Vec(0.3,0.2), poly; include_bounds = false) == true
    # Vec(0.7,0.2)
    @test in(Vec(0.7,0.2), poly; include_bounds = true)  == false
    @test in(Vec(0.7,0.2), poly; include_bounds = false) == false
    # Vec(0.1,0.3)
    @test in(Vec(0.1,0.3), poly; include_bounds = true)  == false
    @test in(Vec(0.1,0.3), poly; include_bounds = false) == false
    # Vec(0.2,0.3)
    @test in(Vec(0.2,0.3), poly; include_bounds = true)  == true
    @test in(Vec(0.2,0.3), poly; include_bounds = false) == true
    # Vec(0.5,0.3)
    @test in(Vec(0.5,0.3), poly; include_bounds = true)  == true
    @test in(Vec(0.5,0.3), poly; include_bounds = false) == true
    # Vec(0.7,0.3)
    @test in(Vec(0.7,0.3), poly; include_bounds = true)  == false
    @test in(Vec(0.7,0.3), poly; include_bounds = false) == false
    # Vec(0.1,0.4)
    @test in(Vec(0.1,0.4), poly; include_bounds = true)  == false
    @test in(Vec(0.1,0.4), poly; include_bounds = false) == false
    # Vec(0.2,0.4)
    @test in(Vec(0.2,0.4), poly; include_bounds = true)  == true
    @test in(Vec(0.2,0.4), poly; include_bounds = false) == false
    # Vec(0.3,0.4)
    @test in(Vec(0.3,0.4), poly; include_bounds = true)  == false
    @test in(Vec(0.3,0.4), poly; include_bounds = false) == false
    # Vec(0.6,0.4)
    @test in(Vec(0.6,0.4), poly; include_bounds = true)  == true
    @test in(Vec(0.6,0.4), poly; include_bounds = false) == false

    # translate!, rotate!
    poly = Polygon([Vec(0.0,0.0), Vec(1.0,0.0), Vec(1.0,1.0), Vec(0.0,1.0)])
    rotate!(poly, Vec(0.0,0.0,1.0) * π/4)
    @test centroid(poly) ≈ [0.5,0.5]
    @test poly ≈ [[0.5, 0.5-1/√2], [0.5+1/√2, 0.5], [0.5, 0.5+1/√2], [0.5-1/√2, 0.5]]
    rotate!(poly, Vec(0.0,0.0,1.0) * π/2)
    @test centroid(poly) ≈ [0.5,0.5]
    @test poly ≈ [[0.5+1/√2, 0.5], [0.5, 0.5+1/√2], [0.5-1/√2, 0.5], [0.5, 0.5-1/√2]]
    xc = centroid(poly)
    for i in eachindex(poly)
        poly[i] = xc + rotate(poly[i] - xc, inv(poly.q))
    end
    @test poly ≈ [[0.0,0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0]]
    translate!(poly, Vec(0.4, 0.3))
    @test poly ≈ [[0.0,0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0]] .+ Vec(0.4, 0.3)

    # distance
    poly = Polygon(Vec{2,Float64}[(0,0), (10,0), (10,10), (0,10)])
    @test distance(poly, Vec(5.0,-1.0), 2) ≈ [0,1]
    @test distance(poly, Vec(11.0,-1.0), 2) ≈ [-1,1]
    @test distance(poly, Vec(11.0,5.0), 2) ≈ [-1,0]
    @test distance(poly, Vec(11.0,11.0), 2) ≈ [-1,-1]
    @test distance(poly, Vec(5.0,11.0), 2) ≈ [0,-1]
    @test distance(poly, Vec(-1.0,11.0), 2) ≈ [1,-1]
    @test distance(poly, Vec(-1.0,5.0), 2) ≈ [1,0]
    @test distance(poly, Vec(-1.0,-1.0), 2) ≈ [1,1]
    # reverse version
    poly = Polygon(Vec{2,Float64}[(0,0), (0,10), (10,10), (10,0)])
    @test distance(poly, Vec(5.0,0.5), 1) ≈ [0,-0.5]
    @test distance(poly, Vec(9.0,0.5), 1) ≈ [0,-0.5]
    @test distance(poly, Vec(9.5,5.0), 1) ≈ [0.5,0]
    @test distance(poly, Vec(9.5,9.0), 1) ≈ [0.5,0]
    @test distance(poly, Vec(5.0,9.5), 1) ≈ [0,0.5]
    @test distance(poly, Vec(1.0,9.5), 1) ≈ [0,0.5]

    # area
    poly = Polygon(Vec{2,Float64}[(0,0), (10,0), (8,5), (3,5)])
    @test area(poly) ≈ (10 + 5) * 5 / 2

    # enlarge
    poly = Polygon(Vec{2,Float64}[(0,0), (1,0), (1,1), (0,1)])
    @test (@inferred enlarge(poly, 1.1))::Polygon ≈ [[-0.05, -0.05], [1.05, -0.05], [1.05, 1.05], [-0.05, 1.05]]
end
