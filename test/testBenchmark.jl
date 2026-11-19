@testset "Benchmarks{$T}" for T in (Float64, Float32)
    suite = vecgeom_suite(T)
    @test length(suite) == 24
    @test all(b -> isfinite(capacity(b.shape)) || b.shape isa MultiUnion, suite)
    # the point sets follow the VecGeom methodology
    for b in suite[[1, 9, 17, end]]
        rng = Geom4hep.Random.MersenneTwister(3)
        tr, n = b.placement, 1000
        pt, dt = Geom4hep._toin_points(rng, b.shape, tr, b.world, n, 0.8)
        @test all(p -> inside(b.shape, tr * p) == kOutside, pt)                # outside the shape
        @test count(i -> distance_to_in(b.shape, tr * pt[i], tr * dt[i]) < Inf, 1:n) == round(Int, 0.8n)
        po, dout = Geom4hep._toout_points(rng, b.shape, n)
        @test all(p -> inside(b.shape, p) == kInside, po)
        @test all(d -> norm(d) ≈ 1, dout)
    end
    rows = benchmark_shape(suite[1]; npoints=64, nrep=1, verbose=false)
    @test length(rows) == 12 && all(r -> r.ns > 0, rows)
    @test Set(r.op for r in rows) == Set((:inside, :contains, :distance_to_in, :safety_to_in, :distance_to_out, :safety_to_out))
end
