using Random, InteractiveUtils

# A representative set of shapes, including phi sections, tilted traps, cut tubes, polycones with
# steps and boolean combinations
function test_shapes(::Type{T}) where T
    V(x, y, z) = Vector3{T}(x, y, z)
    Any[Box{T}(1, 2, 3),
        Tube{T}(1, 2, 3, 0, 2π), Tube{T}(1, 2, 3, 0.3, 2.0),
        Geom4hep.Cone{T}(1, 2, 1, 3, 4, 0, 2π), Geom4hep.Cone{T}(0, 2, 1, 3, 4, 0, π/2),
        Trd{T}(1, 2, 3, 4, 5),
        Trap{T}(2, 0.2, 0.3, 1, 1, 1.2, 0.1, 1.5, 1.4, 1.6, 0.1),
        CutTube{T}(1, 2, 3, 0, 2π, V(0, -0.3, -1), V(0.3, 0, 1)),
        CutTube{T}(1, 2, 3, 0.5, 2.5, V(0, -0.3, -1), V(0.3, 0, 1)),
        Polycone{T}([0., 1, 1, 0.5], [2., 3, 2, 2.5], [-1., 0, 1, 2], 0, 2π),
        Polycone{T}([0., 0, 1, 1, 0.5, 0.5], [2., 2, 3, 3, 1.5, 2.5], [-1., 0, 0, 1, 1, 2], 0.3, 4.0),
        BooleanUnion(Box{T}(1, 1, 1), Tube{T}(0, 1, 2, 0, 2π), Transformation3D{T}(1, 0, 0)),
        BooleanSubtraction(Box{T}(1, 1, 1), Tube{T}(0, 0.5, 2, 0, 2π)),
        BooleanIntersection(Box{T}(1, 1, 1), Geom4hep.Cone{T}(0, 1.5, 0, 0.5, 1.2, 0, 2π)),
        MultiUnion{T}([Box{T}(1, 1, 1), Tube{T}(0, 1, 2, 0, 2π), Box{T}(0.5, 2, 0.5)],
                      [one(Transformation3D{T}), Transformation3D{T}(1.5, 0, 0), Transformation3D{T}(0, 1, 0)])]
end
shapename(s) = string(nameof(typeof(s)))

@testset "Normals{$T}" for T in (Float64, Float32)
    V(x, y, z) = Vector3{T}(x, y, z)
    P(x, y, z) = Point3{T}(x, y, z)
    # Smallest step along the normal that reliably leaves the surface band (± tolerance/2). Points
    # closer than `step` to a sharp edge can legitimately fail the outward/inward check, so a small
    # fraction of failures is accepted.
    step = 2tolerance(T)
    rng = MersenneTwister(2024)

    @testset "outward unit normal: $(shapename(s))" for s in test_shapes(T)
        lo, hi = extent(s)
        nchecked, ngood = 0, 0
        while nchecked < 200
            p = Point3{T}(lo .+ (hi .- lo) .* rand(rng, T, 3))
            inside(s, p) == kInside || continue
            d = V(normalize(randn(rng, 3))...)
            dist = distance_to_out(s, p, d)
            (isfinite(dist) && dist > 0) || continue
            q = p + dist * d
            inside(s, q) == kSurface || continue
            nchecked += 1
            n = normal(s, q)
            @test n isa Vector3{T}
            @test norm(n) ≈ 1 atol = 10eps(T)
            ngood += inside(s, q + step * n) == kOutside && inside(s, q - step * n) != kOutside
        end
        @test ngood >= 0.98 * nchecked
    end

    @testset "not on surface" begin
        @test normal(Box{T}(1, 2, 3), P(0, 0, 0)) == zero(Vector3{T})
        @test normal(Tube{T}(1, 2, 3, 0, 2π), P(0, 0, 0)) == zero(Vector3{T})
        @test normal(Geom4hep.Cone{T}(1, 2, 1, 3, 4, 0, 2π), P(10, 0, 0)) == zero(Vector3{T})
        @test normal(Trap{T}(1, 2, 3, 4, 5), P(0, 0, 0)) == zero(Vector3{T})
    end

    @testset "specific normals" begin
        c = Geom4hep.Cone{T}(1, 2, 1, 3, 4, 0, 2π)
        @test normal(c, P(1.5, 0, 4)) ≈ V(0, 0, 1)
        @test normal(c, P(1.5, 0, -4)) ≈ V(0, 0, -1)
        @test normal(c, P(0, 2.5, 0)) ≈ normalize(V(0, 8, -1))       # outer cone: r = 2.5 + z/8
        @test normal(c, P(-1, 0, 0)) ≈ V(1, 0, 0)                     # inner cylinder (rmin1 = rmin2)
        t = Trap{T}(1, 2, 3, 4, 5)                                     # x: 1 → 2, y: 3 → 4 over z ∈ [-5, 5]
        @test normal(t, P(0, 0, 5)) ≈ V(0, 0, 1)
        @test normal(t, P(1.5, 0, 0)) ≈ normalize(V(10, 0, -1))
        # polycone with a step at z = 0: rmax 2 below, 3 above → ring face looking down
        pc = Polycone{T}([0., 0, 0, 0], [2., 2, 3, 3], [-1., 0, 0, 1], 0, 2π)
        @test normal(pc, P(2.5, 0, 0)) ≈ V(0, 0, -1)
        @test normal(pc, P(0, 0, 1)) ≈ V(0, 0, 1)
        @test normal(pc, P(0, 2, -0.5)) ≈ V(0, 1, 0)
        # subtraction: the surface of the hole looks towards the axis
        s = BooleanSubtraction(Box{T}(1, 1, 1), Tube{T}(0, 0.5, 2, 0, 2π))
        @test normal(s, P(0.5, 0, 0)) ≈ V(-1, 0, 0)
        @test normal(s, P(1, 0, 0)) ≈ V(1, 0, 0)
        u = BooleanUnion(Box{T}(1, 1, 1), Box{T}(1, 1, 1), Transformation3D{T}(1, 0, 0))
        @test normal(u, P(2, 0, 0)) ≈ V(1, 0, 0)
        @test normal(u, P(0.5, 0, 0)) == zero(Vector3{T})           # inside
        @test normal(u, P(0.5, 1, 0)) ≈ V(0, 1, 0)                   # overlapping faces
    end
end

@testset "Capacity and surface" begin
    T = Float64
    @testset "$(shapename(s))" for s in test_shapes(T)
        if s isa Union{Geom4hep.BooleanShape, MultiUnion}
            @test capacity(s; npoints=200_000) ≈ estimate_capacity(s; npoints=200_000)
            @test surface(s; npoints=200_000) > 0
        else
            @test estimate_capacity(s; npoints=200_000) ≈ capacity(s) rtol = 2e-3
            @test estimate_surface(s; npoints=1_000_000) ≈ surface(s) rtol = 2e-2
        end
    end
    # union of two disjoint boxes
    u = BooleanUnion(Box{T}(1, 1, 1), Box{T}(1, 2, 3), Transformation3D{T}(5, 0, 0))
    @test capacity(u) ≈ 8 + 48 rtol = 1e-3
    @test surface(u) ≈ 24 + 88 rtol = 2e-2
    # the extent of a tilted, skewed trap contains all its vertices
    t = Trap{T}(2, 0.2, 0.3, 1, 1, 1.2, 0.1, 1.5, 1.4, 1.6, 0.1)
    lo, hi = extent(t)
    @test all(all(lo .<= v .<= hi) for v in t.coordinates)
end

@testset "Float32 type stability" begin
    T = Float32
    P, V = Point3{T}, Vector3{T}
    funcs = [(inside, (P,)), (safety_to_in, (P,)), (safety_to_out, (P,)), (distance_to_in, (P, V)),
             (distance_to_out, (P, V)), (normal, (P,)), (extent, ()), (capacity, ()), (surface, ())]
    @testset "$(shapename(s))" for s in test_shapes(T)
        for (f, args) in funcs
            s isa Union{Geom4hep.BooleanShape, MultiUnion} && f in (capacity, surface) && continue   # sampling estimates
            sig = (typeof(s), args...)
            rt = only(Base.return_types(f, sig))
            @test isconcretetype(rt)
            @test !occursin("Float64", sprint(io -> code_warntype(io, f, sig; optimize=true)))
        end
    end
end
