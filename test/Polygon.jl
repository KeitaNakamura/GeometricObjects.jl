@testset "Polygon" begin
    poly = Polygon([Vec(0.1,0.1), Vec(0.3,0.1), Vec(0.4,0.0), Vec(0.5,0.1), Vec(0.7,0.1),
                    Vec(0.6,0.4), Vec(0.5,0.4), Vec(0.4,0.3), Vec(0.3,0.3), Vec(0.2,0.4)])
    # Vec(0.0,0.0)
    @test isinside(poly, Vec(0.0,0.0); include_bounds = true)  == false
    @test isinside(poly, Vec(0.0,0.0); include_bounds = false) == false
    # Vec(0.0,0.1)
    @test isinside(poly, Vec(0.0,0.1); include_bounds = true)  == false
    @test isinside(poly, Vec(0.0,0.1); include_bounds = false) == false
    # Vec(0.2,0.1)
    @test isinside(poly, Vec(0.2,0.1); include_bounds = true)  == true
    @test isinside(poly, Vec(0.2,0.1); include_bounds = false) == false
    # Vec(0.4,0.1)
    @test isinside(poly, Vec(0.4,0.1); include_bounds = true)  == true
    @test isinside(poly, Vec(0.4,0.1); include_bounds = false) == true
    # Vec(0.6,0.1)
    @test isinside(poly, Vec(0.6,0.1); include_bounds = true)  == true
    @test isinside(poly, Vec(0.6,0.1); include_bounds = false) == false
    # Vec(0.1,0.2)
    @test isinside(poly, Vec(0.1,0.2); include_bounds = true)  == false
    @test isinside(poly, Vec(0.1,0.2); include_bounds = false) == false
    # Vec(0.3,0.2)
    @test isinside(poly, Vec(0.3,0.2); include_bounds = true)  == true
    @test isinside(poly, Vec(0.3,0.2); include_bounds = false) == true
    # Vec(0.7,0.2)
    @test isinside(poly, Vec(0.7,0.2); include_bounds = true)  == false
    @test isinside(poly, Vec(0.7,0.2); include_bounds = false) == false
    # Vec(0.1,0.3)
    @test isinside(poly, Vec(0.1,0.3); include_bounds = true)  == false
    @test isinside(poly, Vec(0.1,0.3); include_bounds = false) == false
    # Vec(0.2,0.3)
    @test isinside(poly, Vec(0.2,0.3); include_bounds = true)  == true
    @test isinside(poly, Vec(0.2,0.3); include_bounds = false) == true
    # Vec(0.5,0.3)
    @test isinside(poly, Vec(0.5,0.3); include_bounds = true)  == true
    @test isinside(poly, Vec(0.5,0.3); include_bounds = false) == true
    # Vec(0.7,0.3)
    @test isinside(poly, Vec(0.7,0.3); include_bounds = true)  == false
    @test isinside(poly, Vec(0.7,0.3); include_bounds = false) == false
    # Vec(0.1,0.4)
    @test isinside(poly, Vec(0.1,0.4); include_bounds = true)  == false
    @test isinside(poly, Vec(0.1,0.4); include_bounds = false) == false
    # Vec(0.2,0.4)
    @test isinside(poly, Vec(0.2,0.4); include_bounds = true)  == true
    @test isinside(poly, Vec(0.2,0.4); include_bounds = false) == false
    # Vec(0.3,0.4)
    @test isinside(poly, Vec(0.3,0.4); include_bounds = true)  == false
    @test isinside(poly, Vec(0.3,0.4); include_bounds = false) == false
    # Vec(0.6,0.4)
    @test isinside(poly, Vec(0.6,0.4); include_bounds = true)  == true
    @test isinside(poly, Vec(0.6,0.4); include_bounds = false) == false
end
