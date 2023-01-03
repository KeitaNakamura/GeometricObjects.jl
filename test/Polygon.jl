@testset "Polygon" begin
    # `in`
    poly = Polygon(Vec(0.1,0.1), Vec(0.3,0.1), Vec(0.4,0.0), Vec(0.5,0.1), Vec(0.7,0.1),
                   Vec(0.6,0.4), Vec(0.5,0.4), Vec(0.4,0.3), Vec(0.3,0.3), Vec(0.2,0.4))
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

    # findall for `in`
    poly = Polygon(Vec(0.0,0.0), Vec(1.0,0.0), Vec(1.0,1.0), Vec(0.0,1.0))
    points = Vec{2,Float64}[(0.5,0.5), (0.8,0.6), (1.2,1.8), (0.1,0.2)]
    @test findall(in(poly), points) == [1,2,4]

    # rotate
    poly = Polygon(Vec(0.0,0.0), Vec(1.0,0.0), Vec(1.0,1.0), Vec(0.0,1.0))
    poly = @inferred rotate(poly, π/4)
    @test centroid(poly) ≈ [0.5,0.5]
    @test coordinates(poly) ≈ [[0.5, 0.5-1/√2], [0.5+1/√2, 0.5], [0.5, 0.5+1/√2], [0.5-1/√2, 0.5]]
    poly = @inferred rotate(poly, π/2)
    @test centroid(poly) ≈ [0.5,0.5]
    @test coordinates(poly) ≈ [[0.5+1/√2, 0.5], [0.5, 0.5+1/√2], [0.5-1/√2, 0.5], [0.5, 0.5-1/√2]]
    # translate
    poly = Polygon(Vec(0.0,0.0), Vec(1.0,0.0), Vec(1.0,1.0), Vec(0.0,1.0))
    @test coordinates(@inferred translate(poly, Vec(0.4, 0.3))) ≈ coordinates(poly) .+ Vec(0.4, 0.3)

    # distance
    poly = Polygon(Vec{2,Float64}[(0,0), (10,0), (10,10), (0,10)]...)
    @test distance(poly, Vec(5.0,-1.0), 2.0) ≈ [0,1]
    @test distance(poly, Vec(11.0,-1.0), 2.0) ≈ [-1,1]
    @test distance(poly, Vec(11.0,5.0), 2.0) ≈ [-1,0]
    @test distance(poly, Vec(11.0,11.0), 2.0) ≈ [-1,-1]
    @test distance(poly, Vec(5.0,11.0), 2.0) ≈ [0,-1]
    @test distance(poly, Vec(-1.0,11.0), 2.0) ≈ [1,-1]
    @test distance(poly, Vec(-1.0,5.0), 2.0) ≈ [1,0]
    @test distance(poly, Vec(-1.0,-1.0), 2.0) ≈ [1,1]
    @test distance(poly, Vec(5.0,-1.0), 2.0, :return_indices) == ([0,1], [1])
    @test distance(poly, Vec(11.0,-1.0), 2.0, :return_indices) == ([-1,1], [1,2])
    @test distance(poly, Vec(11.0,5.0), 2.0, :return_indices) == ([-1,0], [2])
    @test distance(poly, Vec(11.0,11.0), 2.0, :return_indices) == ([-1,-1], [2,3])
    @test distance(poly, Vec(5.0,11.0), 2.0, :return_indices) == ([0,-1], [3])
    @test distance(poly, Vec(-1.0,11.0), 2.0, :return_indices) == ([1,-1], [3,4])
    @test distance(poly, Vec(-1.0,5.0), 2.0, :return_indices) == ([1,0], [4])
    @test distance(poly, Vec(-1.0,-1.0), 2.0, :return_indices) == ([1,1], [4,1])
    # reverse version
    poly = Polygon(Vec{2,Float64}[(0,0), (0,10), (10,10), (10,0)]...)
    @test distance(poly, Vec(5.0,0.5), 1.0) ≈ [0,-0.5]
    @test distance(poly, Vec(9.0,0.5), 1.0) ≈ [0.5,-0.25]
    @test distance(poly, Vec(9.5,5.0), 1.0) ≈ [0.5,0]
    @test distance(poly, Vec(9.5,9.0), 1.0) ≈ [0.25,0.5]
    @test distance(poly, Vec(5.0,9.5), 1.0) ≈ [0,0.5]
    @test distance(poly, Vec(1.0,9.5), 1.0) ≈ [-0.5,0.25]
    @test distance(poly, Vec(5.0,0.5), 1.0, :return_indices) == ([0,-0.5], [4])
    @test distance(poly, Vec(9.0,0.5), 1.0, :return_indices) == ([0.5,-0.25], [3,4])
    @test distance(poly, Vec(9.5,5.0), 1.0, :return_indices) == ([0.5,0], [3])
    @test distance(poly, Vec(9.5,9.0), 1.0, :return_indices) == ([0.25,0.5], [2,3])
    @test distance(poly, Vec(5.0,9.5), 1.0, :return_indices) == ([0,0.5], [2])
    @test distance(poly, Vec(1.0,9.5), 1.0, :return_indices) == ([-0.5,0.25], [1,2])

    # area
    poly = Polygon(Vec{2,Float64}[(0,0), (10,0), (8,5), (3,5)]...)
    @test (@inferred area(poly))::Float64 ≈ (10 + 5) * 5 / 2

    # enlarge
    poly = Polygon(Vec{2,Float64}[(0,0), (1,0), (1,1), (0,1)]...)
    @test (@inferred enlarge(poly, 1.1))::Polygon |> coordinates ≈ [[-0.05, -0.05], [1.05, -0.05], [1.05, 1.05], [-0.05, 1.05]]

    # centered
    poly = Polygon(Vec{2,Float64}[(0,0), (1,0), (1,1), (0,1)]...)
    @test (@inferred GeometricObjects.centered(poly))::Polygon |> coordinates ≈ [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]

    # intersect
    poly = Polygon(Vec(0.0,1.0), Vec(1.0,0.0), Vec(1.0,1.0))
    line = Line(Vec(0.0,0.0), Vec(0.2,0.2));
    @test intersect(poly, line, extended = true) ≈ [0.5,0.5]
    @test intersect(poly, line, extended = false) === nothing
    line = Line(Vec(2.0,-0.2), Vec(1.8,-0.2));
    @test intersect(poly, line, extended = true) === nothing
    @test intersect(poly, line, extended = false) === nothing

    # moment of inertia
    poly = Rectangle(Vec(-2.0, -1.0), Vec(1.0, 3.0))
    @test (@inferred moment_of_inertia(poly)) ≈ (3^2+4^2)/12
end
