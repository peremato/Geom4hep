using Random

@testset "Torus{$T}" for T in (Float64, Float32)
    tol = Geom4hep.tolerance(T)
    vx, vy, vz = Vector3{T}(1, 0, 0), Vector3{T}(0, 1, 0), Vector3{T}(0, 0, 1)
    fat    = Torus{T}(0, 20, 50)
    hollow = Torus{T}(10, 20, 50)
    seg    = Torus{T}(0, 20, 50, 0, π/2)                     # Δϕ ≤ π
    big    = Torus{T}(5, 20, 50, π/4, 3π/2)                  # Δϕ > π, hollow
    thin   = Torus{T}(0, 12.5, 1326)                         # as in CMS (rtor/rmax ~ 100)
    thinseg = Torus{T}(0, 7, 1674.2, deg2rad(2.5), deg2rad(85))
    shapes = (fat, hollow, seg, big, thin, thinseg)

    @testset "construction" begin
        @test isbitstype(typeof(seg))
        @test_throws ErrorException Torus{T}(0, 60, 50)       # rmax ≥ rtor
        @test_throws ErrorException Torus{T}(20, 10, 50)
        @test all(extent(fat) .≈ (Point3{T}(-70, -70, -20), Point3{T}(70, 70, 20)))
        lo, hi = extent(seg)
        @test isapprox(lo, Point3{T}(0, 0, -20); atol=1e-4) && hi ≈ Point3{T}(70, 70, 20)
    end

    @testset "inside" begin
        @test inside(fat, Point3{T}(50, 0, 0)) == kInside
        @test inside(fat, Point3{T}(0, 0, 0)) == kOutside       # the hole
        @test inside(fat, Point3{T}(70, 0, 0)) == kSurface
        @test inside(fat, Point3{T}(0, 50, 20)) == kSurface
        @test inside(fat, Point3{T}(0, 50, 21)) == kOutside
        @test inside(hollow, Point3{T}(50, 0, 0)) == kOutside   # the core of the tube
        @test inside(hollow, Point3{T}(65, 0, 0)) == kInside
        @test inside(hollow, Point3{T}(60, 0, 0)) == kSurface
        @test inside(seg, Point3{T}(0, -50, 0)) == kOutside
        @test inside(seg, Point3{T}(50, 0, 0)) == kSurface       # on the ϕ = 0 plane
        @test inside(big, Point3{T}(-65, 0, 0)) == kInside
        @test inside(big, Point3{T}(65, 1, 0)) == kOutside       # in the missing wedge
        @test inside(thin, Point3{T}(0, -1326, 12)) == kInside
        @test inside(thin, Point3{T}(0, -1326, 13)) == kOutside
    end

    @testset "distances" begin
        @test distance_to_in(fat, Point3{T}(-100, 0, 0), vx) ≈ 30
        @test distance_to_in(fat, Point3{T}(0, 0, 0), vx) ≈ 30                 # from the hole
        @test distance_to_in(fat, Point3{T}(0, 0, -100), vz) == Inf            # through the hole
        @test distance_to_in(fat, Point3{T}(50, 0, -100), vz) ≈ 80
        @test distance_to_in(fat, Point3{T}(-100, 69, 0), vx) ≈ 100 - sqrt(T(70^2 - 69^2)) rtol=1e-4   # grazing
        @test distance_to_in(fat, Point3{T}(-100, 71, 0), vx) == Inf
        @test distance_to_out(fat, Point3{T}(50, 0, 0), vx) ≈ 20
        @test distance_to_out(fat, Point3{T}(50, 0, 0), vz) ≈ 20
        @test distance_to_out(fat, Point3{T}(-50, 0, 0), vx) ≈ 20               # into the hole
        @test distance_to_in(hollow, Point3{T}(50, 0, 0), vx) ≈ 10              # from the core
        @test distance_to_out(hollow, Point3{T}(65, 0, 0), -vx) ≈ 5             # into the core
        @test distance_to_in(seg, Point3{T}(50, -10, 0), vy) ≈ 10               # through ϕ = 0
        @test distance_to_out(seg, Point3{T}(50, 10, 0), -vy) ≈ 10
        @test distance_to_in(big, Point3{T}(40, 1, 0), vy) ≈ 39 rtol=1e-4      # out of the missing wedge
        @test distance_to_in(thin, Point3{T}(0, -2000, 1), vy) ≈ 2000 - 1326 - sqrt(T(12.5^2 - 1))
        @test distance_to_in(thin, Point3{T}(-2000, 0, 12.4), vx) ≈ 2000 - 1326 - sqrt(T(12.5^2 - 12.4^2)) rtol=1e-5
        @test distance_to_in(fat, Point3{T}(70, 0, 0), -vx) <= 0                # on the surface, moving in
        @test distance_to_in(fat, Point3{T}(70, 0, 0), vx) == Inf               # moving out
        @test distance_to_out(fat, Point3{T}(70, 0, 0), vx) ≈ 0 atol=tol
        @test distance_to_out(fat, Point3{T}(100, 0, 0), vx) == -1
        @test distance_to_in(fat, Point3{T}(50, 0, 0), vx) == -1
    end

    @testset "normals" begin
        @test normal(fat, Point3{T}(70, 0, 0)) ≈ vx
        @test normal(fat, Point3{T}(30, 0, 0)) ≈ -vx
        @test normal(fat, Point3{T}(0, 50, 20)) ≈ vz
        @test normal(hollow, Point3{T}(60, 0, 0)) ≈ -vx
        @test normal(seg, Point3{T}(50, 0, 0)) ≈ -vy
        @test normal(fat, Point3{T}(50, 0, 0)) == zero(Vector3{T})
    end

    @testset "capacity and surface" begin
        @test capacity(fat) ≈ T(π * 400 * 50 * 2π)
        @test surface(fat) ≈ T(2π * 20 * 50 * 2π)
        @test surface(seg) ≈ T(2π * 20 * 50 * π/2 + 2π * 400)
        for s in (fat, hollow, seg, big)
            @test estimate_capacity(s; npoints=200_000) ≈ capacity(s) rtol=0.03
            @test estimate_surface(s; npoints=400_000) ≈ surface(s) rtol=0.05
        end
    end

    # Random rays compared with sampling of `inside` along the ray: a miss must not cross the shape,
    # an exit must not leave it before, and the end points must be on the surface. For the thin tori
    # the points are taken close to the tube.
    @testset "random rays $(s)" for s in shapes
        rng = MersenneTwister(2468)
        lo, hi = extent(s)
        thintube = s.rmax < s.rtor / 50
        δ = thintube ? s.rmax / 25 : s.rmax / 100
        L = sqrt((hi - lo) ⋅ (hi - lo))
        nin = nout = bad = nhit = 0
        for _ in 1:(thintube ? 600 : 1000)
            if thintube
                ϕ = s.ϕ₀ + s.Δϕ * rand(rng, T)
                p = Point3{T}(s.rtor * cos(ϕ), s.rtor * sin(ϕ), 0) + Vector3{T}((rand(rng, T, 3) .* 2 .- 1) .* 3s.rmax)
            else
                p = Point3{T}((lo + hi) / 2 + (rand(rng, T, 3) .* 2 .- 1) .* 0.6 .* (hi - lo))
            end
            d = normalize(Vector3{T}(randn(rng, T, 3)))
            pos = inside(s, p)
            if pos == kOutside
                nin += 1
                t = distance_to_in(s, p, d)
                if t == Inf
                    bad += any(k -> inside(s, p + (k * δ) * d) == kInside, 1:ceil(Int, L / δ))
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
        @test nhit > 30 && nout > 30
        @test bad <= 0.005 * (nin + nout)          # a few rays may graze the surfaces
    end

    @testset "mesh" begin
        function meshvolume(s, facets)
            o = Vector3{Float64}(1, 2, 3)
            p = collect(coordinates(s, facets))
            sum(faces(s, facets)) do f
                a, b, c, d = (Vector3{Float64}(p[k]) - o for k in f)
                (a ⋅ (b × c) + a ⋅ (c × d)) / 6
            end
        end
        for s in (fat, hollow, seg, big, thinseg)
            @test meshvolume(s, 400) ≈ capacity(s) rtol = 2e-3
        end
    end
end
