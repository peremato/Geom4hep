@testset "CutTube{$T}" for T in (Float64, Float32)
    P(x, y, z) = Point3{T}(x, y, z)
    V(x, y, z) = Vector3{T}(x, y, z)

    # With planes perpendicular to z the cut tube is an ordinary tube
    ct = CutTube{T}(5, 10, 20, 0, 2π, V(0, 0, -1), V(0, 0, 1))
    tb = Tube{T}(5, 10, 20, 0, 2π)
    for p in (P(7, 0, 0), P(0, 0, 0), P(7, 0, 19), P(7, 0, 25), P(0, 12, 0), P(10, 0, 0), P(7, 0, 20))
        @test inside(ct, p) == inside(tb, p)
    end
    @test distanceToOut(ct, P(7, 0, 0), V(0, 0, 1)) ≈ 20
    @test distanceToOut(ct, P(7, 0, 0), V(1, 0, 0)) ≈ 3
    @test distanceToIn(ct, P(7, 0, 50), V(0, 0, -1)) ≈ 30
    @test distanceToIn(ct, P(-30, 0, 0), V(1, 0, 0)) ≈ 20
    lo, hi = extent(ct)
    @test lo[3] ≈ -20 && hi[3] ≈ 20
    @test hi[1] ≈ 10 && hi[2] ≈ 10

    # Tilted top plane (45° around y): top surface is z = 20 - x
    s = 1/sqrt(T(2))
    ct45 = CutTube{T}(0, 10, 20, 0, 2π, V(0, 0, -1), V(s, 0, s))
    @test inside(ct45, P(5, 0, 14)) == kInside
    @test inside(ct45, P(5, 0, 16)) == kOutside
    @test inside(ct45, P(-5, 0, 24)) == kInside
    @test distanceToOut(ct45, P(5, 0, 0), V(0, 0, 1)) ≈ 15
    @test distanceToOut(ct45, P(-5, 0, 0), V(0, 0, 1)) ≈ 25
end
