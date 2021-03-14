@testset "Polygon" begin
    poly = Polygon([Vec(0.1,0.1), Vec(0.3,0.1), Vec(0.4,0.0), Vec(0.5,0.1), Vec(0.7,0.1),
                    Vec(0.6,0.4), Vec(0.5,0.4), Vec(0.4,0.3), Vec(0.3,0.3), Vec(0.2,0.4)])
    # Vec(0.0,0.0)
    @test within(Vec(0.0,0.0), poly; include_bounds = true)  == false
    @test within(Vec(0.0,0.0), poly; include_bounds = false) == false
    # Vec(0.0,0.1)
    @test within(Vec(0.0,0.1), poly; include_bounds = true)  == false
    @test within(Vec(0.0,0.1), poly; include_bounds = false) == false
    # Vec(0.2,0.1)
    @test within(Vec(0.2,0.1), poly; include_bounds = true)  == true
    @test within(Vec(0.2,0.1), poly; include_bounds = false) == false
    # Vec(0.4,0.1)
    @test within(Vec(0.4,0.1), poly; include_bounds = true)  == true
    @test within(Vec(0.4,0.1), poly; include_bounds = false) == true
    # Vec(0.6,0.1)
    @test within(Vec(0.6,0.1), poly; include_bounds = true)  == true
    @test within(Vec(0.6,0.1), poly; include_bounds = false) == false
    # Vec(0.1,0.2)
    @test within(Vec(0.1,0.2), poly; include_bounds = true)  == false
    @test within(Vec(0.1,0.2), poly; include_bounds = false) == false
    # Vec(0.3,0.2)
    @test within(Vec(0.3,0.2), poly; include_bounds = true)  == true
    @test within(Vec(0.3,0.2), poly; include_bounds = false) == true
    # Vec(0.7,0.2)
    @test within(Vec(0.7,0.2), poly; include_bounds = true)  == false
    @test within(Vec(0.7,0.2), poly; include_bounds = false) == false
    # Vec(0.1,0.3)
    @test within(Vec(0.1,0.3), poly; include_bounds = true)  == false
    @test within(Vec(0.1,0.3), poly; include_bounds = false) == false
    # Vec(0.2,0.3)
    @test within(Vec(0.2,0.3), poly; include_bounds = true)  == true
    @test within(Vec(0.2,0.3), poly; include_bounds = false) == true
    # Vec(0.5,0.3)
    @test within(Vec(0.5,0.3), poly; include_bounds = true)  == true
    @test within(Vec(0.5,0.3), poly; include_bounds = false) == true
    # Vec(0.7,0.3)
    @test within(Vec(0.7,0.3), poly; include_bounds = true)  == false
    @test within(Vec(0.7,0.3), poly; include_bounds = false) == false
    # Vec(0.1,0.4)
    @test within(Vec(0.1,0.4), poly; include_bounds = true)  == false
    @test within(Vec(0.1,0.4), poly; include_bounds = false) == false
    # Vec(0.2,0.4)
    @test within(Vec(0.2,0.4), poly; include_bounds = true)  == true
    @test within(Vec(0.2,0.4), poly; include_bounds = false) == false
    # Vec(0.3,0.4)
    @test within(Vec(0.3,0.4), poly; include_bounds = true)  == false
    @test within(Vec(0.3,0.4), poly; include_bounds = false) == false
    # Vec(0.6,0.4)
    @test within(Vec(0.6,0.4), poly; include_bounds = true)  == true
    @test within(Vec(0.6,0.4), poly; include_bounds = false) == false
end
