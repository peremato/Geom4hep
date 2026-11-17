using Random

@testset "Polyhedra{$T}" for T in (Float64, Float32)
    tol = Geom4hep.tolerance(T)
    vx, vy, vz = Vector3{T}(1, 0, 0), Vector3{T}(0, 1, 0), Vector3{T}(0, 0, 1)
    c30, s30 = cos(T(π)/6), sin(T(π)/6)
    hexa  = Polyhedra{T}(6, [0, 0], [10, 10], [-5, 5])                   # hexagonal prism, corners on x
    hollow = Polyhedra{T}(8, [5, 5, 8, 8], [20, 20, 20, 15], [-30, 0, 0, 30])  # a step at z = 0
    seg   = Polyhedra{T}(3, [10, 10], [40, 40], [-10, 10], 0, π/2)        # Δϕ ≤ π
    big   = Polyhedra{T}(6, [10, 15, 10], [40, 45, 40], [-20, 0, 20], π/4, 3π/2)   # Δϕ > π, kinks in z
    wedge = Polyhedra{T}(1, [1775, 1775, 1775, 2770.71, 2813.42, 2813.42],      # as in the CMS HCAL
                            [1866.5, 1866.5, 1927.03, 2870, 2870, 2870],
                            [3242, 3704.6, 3750.22, 4460.79, 4491.27, 5541], deg2rad(352.5), deg2rad(15))
    cone  = Polyhedra{T}(5, [0, 0], [0, 30], [0, 30])                     # pentagonal pyramid
    shapes = (hexa, hollow, seg, big, wedge, cone)

    @testset "construction" begin
        @test isbitstype(typeof(hexa))
        @test nsections(hollow) == 2                                      # the zero-length step is dropped
        @test_throws ErrorException Polyhedra{T}(1, [0, 0], [10, 10], [0, 1])         # one side, full ϕ
        @test_throws ErrorException Polyhedra{T}(6, [0, 0, 0], [10, 10, 10], [0, 2, 1]) # z not ordered
        @test all(extent(hexa) .≈ (Point3{T}(-10/c30, -10, -5), Point3{T}(10/c30, 10, 5)))
        lo, hi = extent(seg)
        @test isapprox(lo, Point3{T}(0, 0, -10); atol=1e-4) && hi ≈ Point3{T}(40/cos(T(π)/12), 40/cos(T(π)/12), 10)
        r = Polyhedra{T}(6, [0, 0], [10, 10], [5, -5])                     # z planes in decreasing order
        @test extent(r) == extent(hexa)
    end

    @testset "inside" begin
        @test inside(hexa, Point3{T}(0, 0, 0)) == kInside
        @test inside(hexa, Point3{T}(11, 0, 0)) == kInside                # towards a corner
        @test inside(hexa, Point3{T}(11c30, 11s30, 0)) == kOutside        # beyond a side
        @test inside(hexa, Point3{T}(10c30, 10s30, 0)) == kSurface        # on a side
        @test inside(hexa, Point3{T}(0, 0, 5)) == kSurface
        @test inside(hollow, Point3{T}(0, 0, 0)) == kOutside
        @test inside(hollow, Point3{T}(0, 6, -10)) == kInside
        @test inside(hollow, Point3{T}(0, 6, 10)) == kOutside             # in the larger hole above z = 0
        @test inside(hollow, Point3{T}(0, 6, 0)) == kSurface              # on the step
        @test inside(hollow, Point3{T}(0, 10, 0)) == kInside              # across the shared plane
        @test inside(seg, Point3{T}(20, 20, 0)) == kInside
        @test inside(seg, Point3{T}(-20, 20, 0)) == kOutside
        @test inside(seg, Point3{T}(20, 0, 0)) == kSurface                # on the ϕ = 0 plane
        @test inside(big, Point3{T}(-20, 0, 0)) == kInside                # in the concave ϕ range
        @test inside(big, Point3{T}(20, 1, 0)) == kOutside                # in the missing wedge
        @test inside(cone, Point3{T}(0, 0, 29)) == kInside
        @test inside(cone, Point3{T}(10, 0, 5)) == kOutside
    end

    @testset "distances" begin
        @test distance_to_in(hexa, Point3{T}(-100, 0, 0), vx) ≈ 100 - 10/c30          # through a corner
        @test distance_to_in(hexa, Point3{T}(0, -100, 0), vy) ≈ 90                     # through a side
        @test distance_to_out(hexa, Point3{T}(0, 0, 0), vy) ≈ 10
        @test distance_to_out(hexa, Point3{T}(0, 0, 0), vz) ≈ 5
        @test distance_to_in(hollow, Point3{T}(0, 0, -10), Vector3{T}(cos(π/8), sin(π/8), 0)) ≈ 5   # from the hole
        @test distance_to_out(hollow, Point3{T}(0, 6, -10), vz) ≈ 10                   # up to the step
        @test distance_to_out(hollow, Point3{T}(0, 10, -10), vz) ≈ 40                  # across z = 0
        @test distance_to_in(hollow, Point3{T}(0, 6, 10), -vz) ≈ 10                    # onto the step
        @test distance_to_in(seg, Point3{T}(20, -20, 0), vy) ≈ 20                      # through ϕ = 0
        @test distance_to_in(big, Point3{T}(40, 30, 0), -vx) ≈ 10                      # from the missing wedge
        @test distance_to_in(hexa, Point3{T}(10c30, 10s30, 0), -Vector3{T}(c30, s30, 0)) <= 0    # on the surface, moving in
        @test distance_to_in(hexa, Point3{T}(10c30, 10s30, 0), Vector3{T}(c30, s30, 0)) == Inf  # moving out
        @test distance_to_out(hexa, Point3{T}(10c30, 10s30, 0), Vector3{T}(c30, s30, 0)) ≈ 0 atol=tol
        @test distance_to_out(hexa, Point3{T}(100, 0, 0), vx) == -1
        @test distance_to_in(hexa, Point3{T}(0, 0, 0), vx) == -1
    end

    @testset "normals" begin
        @test normal(hexa, Point3{T}(10c30, 10s30, 0)) ≈ Vector3{T}(c30, s30, 0)
        @test normal(hexa, Point3{T}(0, 0, 5)) ≈ vz
        @test normal(hollow, Point3{T}(0, 6, 0)) ≈ vz                     # the step faces up
        @test normal(hollow, Point3{T}(5cos(π/8), 5sin(π/8), -10)) ≈ -Vector3{T}(cos(π/8), sin(π/8), 0)
        @test normal(seg, Point3{T}(20, 0, 0)) ≈ -vy
        @test normal(hexa, Point3{T}(0, 0, 0)) == zero(Vector3{T})
        @test normal(hollow, Point3{T}(0, 10, 0)) == zero(Vector3{T})    # shared plane, not a boundary
    end

    @testset "capacity and surface" begin
        @test capacity(hexa) ≈ 6 * tan(T(π)/6) * 100 * 10
        @test capacity(cone) ≈ 5 * tan(T(π)/5) * 30^2 * 30 / 3
        for s in shapes
            @test estimate_capacity(s; npoints=200_000) ≈ capacity(s) rtol=0.03
        end
        @test surface(hexa) ≈ 2 * 6 * tan(T(π)/6) * 100 + 6 * 2 * 10 * tan(T(π)/6) * 10
        for s in (hexa, hollow, seg, big, cone)
            @test estimate_surface(s; npoints=400_000) ≈ surface(s) rtol=0.05
        end
    end

    # Random rays compared with sampling of `inside` along the ray: a miss must not cross the shape,
    # an exit must not leave it before, and the end points must be on the surface
    @testset "random rays $(s)" for s in shapes
        rng = MersenneTwister(4321)
        lo, hi = extent(s)
        c, h = (lo + hi) / 2, (hi - lo) / 2
        L = 2 * sqrt(h ⋅ h)
        δ = L / 2000
        nin = nout = bad = nhit = 0
        for _ in 1:(s === wedge ? 5000 : 1000)          # the wedge fills a small part of its extent
            p = Point3{T}(c + (rand(rng, T, 3) .* 2 .- 1) .* 1.2h)
            d = normalize(Vector3{T}(randn(rng, T, 3)))
            pos = inside(s, p)
            if pos == kOutside
                nin += 1
                t = distance_to_in(s, p, d)
                if t == Inf
                    bad += any(k -> inside(s, p + (k * δ) * d) == kInside, 1:4000)
                else
                    nhit += 1
                    q = p + t * d
                    bad += !(inside(s, q) == kSurface && inside(s, q + 10tol * d) != kOutside)
                    bad += any(k -> inside(s, p + (k * δ) * d) == kInside, 1:floor(Int, t / δ) - 1)
                end
            elseif pos == kInside
                nout += 1
                t = distance_to_out(s, p, d)
                q = p + t * d
                bad += !(inside(s, q) == kSurface && inside(s, q + 10tol * d) != kInside)
                bad += any(k -> inside(s, p + (k * δ) * d) == kOutside, 1:floor(Int, t / δ) - 1)
            end
        end
        @test nhit > 50 && nout > 50
        @test bad <= 0.005 * (nin + nout)          # a few rays may graze edges
    end

    @testset "mesh" begin
        function meshvolume(s)
            o = Vector3{Float64}(1, 2, 3)
            p = collect(coordinates(s))
            sum(faces(s)) do f
                a, b, c, d = (Vector3{Float64}(p[k]) - o for k in f)
                (a ⋅ (b × c) + a ⋅ (c × d)) / 6
            end
        end
        for s in shapes
            @test meshvolume(s) ≈ capacity(s) rtol = 1e-4
        end
        @test length(coordinates(mesh(hollow))) == 2 * 4 * 9
    end
end
