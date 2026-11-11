using Random

@testset "MultiUnion{$T}" for T in (Float64, Float32)
    P(x, y, z) = Point3{T}(x, y, z)
    V(x, y, z) = Vector3{T}(x, y, z)
    tr(x, y, z) = Transformation3D{T}(x, y, z)
    shapes = [Box{T}(1, 1, 1), Tube{T}(0, 1, 2, 0, 2π), Box{T}(0.5, 2, 0.5), Geom4hep.Cone{T}(0, 1, 0, 0.5, 1, 0, 2π)]
    places = [tr(0, 0, 0), tr(1.5, 0, 0), tr(0, 1, 0), Transformation3D{T}(-1, 0, 0.5, 0, π/2, 0)]

    mu = MultiUnion{T}(shapes, places)
    @test ncomponents(mu) == 4
    # the same shape as a chain of binary unions
    chain = BooleanUnion(shapes[1], shapes[2], places[2])
    chain = BooleanUnion(chain, shapes[3], places[3])
    chain = BooleanUnion(chain, shapes[4], places[4])

    lo, hi = extent(mu)
    clo, chi = extent(chain)
    @test lo ≈ clo && hi ≈ chi

    @test inside(mu, P(0, 0, 0)) == kInside
    @test inside(mu, P(1.5, 0, 1.9)) == kInside
    @test inside(mu, P(0, 3.5, 0)) == kOutside
    @test inside(mu, P(0, 3, 0)) == kSurface
    @test distance_to_in(mu, P(0, 10, 0), V(0, -1, 0)) ≈ 7
    @test distance_to_out(mu, P(0, 0, 0), V(1, 0, 0)) ≈ 2.5 atol = 10push_tolerance(T)  # box, then tube
    @test distance_to_out(mu, P(0, 0, 0), V(0, 1, 0)) ≈ 3 atol = 10push_tolerance(T)
    @test distance_to_out(mu, P(5, 5, 5), V(0, 1, 0)) == -1

    # agreement with the chain of binary unions on random points and rays
    rng = MersenneTwister(11)
    nin, nout = 0, 0
    for _ in 1:2000
        p = Point3{T}(lo .- 1 .+ (hi .- lo .+ 2) .* rand(rng, T, 3))
        d = V(normalize(randn(rng, 3))...)
        pos = inside(mu, p)
        @test pos == inside(chain, p)
        if pos == kInside
            nin += 1
            @test distance_to_out(mu, p, d) ≈ distance_to_out(chain, p, d) atol = 100push_tolerance(T)
        elseif pos == kOutside
            nout += 1
            @test distance_to_in(mu, p, d) ≈ distance_to_in(chain, p, d)
        end
    end
    @test nin > 100 && nout > 100

    @testset "add_component!" begin
        m = MultiUnion{T}()
        add_component!(m, Box{T}(1, 1, 1))
        add_component!(m, Box{T}(1, 1, 1), tr(3, 0, 0))
        @test ncomponents(m) == 2
        @test inside(m, P(3, 0, 0)) == kInside
        @test inside(m, P(1.5, 0, 0)) == kOutside
        @test capacity(m) ≈ 16 rtol = 1e-3
    end
end

@testset "Boolean fixes{$T}" for T in (Float64, Float32)
    P(x, y, z) = Point3{T}(x, y, z)
    V(x, y, z) = Vector3{T}(x, y, z)
    # subtraction with a displaced subtrahend: box [-10,10]^3 minus box [-2,2]^3 shifted to x=+20
    # (outside the box) and a hole box shifted to x = 5
    s = BooleanSubtraction(Box{T}(10, 10, 10), Box{T}(2, 2, 20), Transformation3D{T}(5, 0, 0))
    @test inside(s, P(5, 0, 0)) == kOutside
    @test inside(s, P(0, 0, 0)) == kInside
    @test distance_to_in(s, P(5, 0, 30), V(0, 0, -1)) == T(Inf)           # through the hole
    @test distance_to_in(s, P(0, 0, 30), V(0, 0, -1)) ≈ 20
    @test distance_to_in(s, P(-30, 0, 0), V(1, 0, 0)) ≈ 20
    @test distance_to_in(s, P(5, 0, 0), V(1, 0, 0)) ≈ 2                    # from inside the hole
    # intersection that needs several propagations: two long boxes crossing at an angle
    i = BooleanIntersection(Box{T}(10, 1, 1), Box{T}(1, 10, 1), Transformation3D{T}(3, 0, 0))
    @test inside(i, P(3, 0, 0)) == kInside
    @test inside(i, P(0, 0, 0)) == kOutside
    @test distance_to_in(i, P(-20, 0, 0), V(1, 0, 0)) ≈ 22                # enters left at -10, right at 2
    @test distance_to_in(i, P(3, -20, 0), V(0, 1, 0)) ≈ 19                # enters right at -10, left at -1
end
