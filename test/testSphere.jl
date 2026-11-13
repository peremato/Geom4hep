using Random

@testset "Sphere{$T}" for T in (Float64, Float32)
    tol = Geom4hep.tolerance(T)
    vx, vy, vz = Vector3{T}(1, 0, 0), Vector3{T}(0, 1, 0), Vector3{T}(0, 0, 1)
    orb   = Sphere{T}(0, 50)
    shell = Sphere{T}(40, 50)
    upper = Sphere{T}(40, 50, 0, 2π, 0, π/2)
    lower = Sphere{T}(0, 50, 0, 2π, π/2, π/2)
    seg   = Sphere{T}(40, 50, 0, π/2)
    seg2  = Sphere{T}(0, 50, π/4, 3π/2, 0, π/2)      # concave ϕ wedge and upper hemisphere
    cap   = Sphere{T}(0, 50, 0, 2π, 0, π/6)          # cone θ ≤ 30°
    band  = Sphere{T}(20, 50, 0, 2π, π/4, π/2)       # two cones, 45° ≤ θ ≤ 135°
    low   = Sphere{T}(40, 50, π/4, π, 2π/3, π/4)     # two cones on the lower nappe, ϕ segment
    wide  = Sphere{T}(0, 50, 0, 2π, 0, 2π/3)         # θ ≤ 120°: concave
    upcone = Sphere{T}(0, 50, 0, 3π/2, π/6, π/3)     # cone θ ≥ 30°, z ≥ 0 plane, concave ϕ wedge

    @testset "construction" begin
        @test_throws ErrorException Sphere{T}(0, 50, 0, 2π, π/2, 2π/3)    # θ beyond π
        @test_throws ErrorException Sphere{T}(50, 40)
        @test extent(upper) == (Point3{T}(-50, -50, 0), Point3{T}(50, 50, 50))
        @test extent(lower) == (Point3{T}(-50, -50, -50), Point3{T}(50, 50, 0))
        @test isbitstype(typeof(seg2))
        @test all(extent(cap) .≈ (Point3{T}(-25, -25, 0), Point3{T}(25, 25, 50)))
        @test all(extent(band) .≈ (Point3{T}(-50, -50, -50cos(π/4)), Point3{T}(50, 50, 50cos(π/4))))
        @test all(extent(low) .≈ (Point3{T}(-50sin(2π/3), -50sin(2π/3), 50cos(11π/12)), Point3{T}(50sin(2π/3), 50sin(2π/3), -20)))
    end

    @testset "inside" begin
        @test inside(orb, Point3{T}(0, 0, 0)) == kInside
        @test inside(orb, Point3{T}(50, 0, 0)) == kSurface
        @test inside(orb, Point3{T}(51, 0, 0)) == kOutside
        @test inside(shell, Point3{T}(0, 0, 0)) == kOutside
        @test inside(shell, Point3{T}(0, 45, 0)) == kInside
        @test inside(shell, Point3{T}(0, 40, 0)) == kSurface
        @test inside(upper, Point3{T}(0, 0, 45)) == kInside
        @test inside(upper, Point3{T}(0, 0, -45)) == kOutside
        @test inside(upper, Point3{T}(45, 0, 0)) == kSurface
        @test inside(seg, Point3{T}(30, 30, 0)) == kInside
        @test inside(seg, Point3{T}(-30, 30, 0)) == kOutside
        @test inside(seg, Point3{T}(45, 0, 0)) == kSurface
        @test inside(seg2, Point3{T}(-30, 0, 10)) == kInside      # in the concave wedge
        @test inside(seg2, Point3{T}(30, 1, 10)) == kOutside
        @test inside(cap, Point3{T}(0, 0, 40)) == kInside
        @test inside(cap, Point3{T}(30, 0, 30)) == kOutside
        @test inside(cap, Point3{T}(20, 0, 40cos(π/6))) == kSurface
        @test inside(cap, Point3{T}(0, 0, -10)) == kOutside
        @test inside(band, Point3{T}(30, 0, 10)) == kInside
        @test inside(band, Point3{T}(10, 0, 30)) == kOutside
        @test inside(band, Point3{T}(10, 0, -30)) == kOutside
        @test inside(low, Point3{T}(0, 45sin(7π/9), 45cos(7π/9))) == kInside
        @test inside(low, Point3{T}(0, -45sin(7π/9), 45cos(7π/9))) == kOutside
        @test inside(wide, Point3{T}(10, 0, -5)) == kInside
        @test inside(wide, Point3{T}(1, 0, -30)) == kOutside
    end

    @testset "distances" begin
        @test distance_to_in(orb, Point3{T}(-100, 0, 0), vx) ≈ 50
        @test distance_to_in(orb, Point3{T}(-100, 60, 0), vx) == Inf
        @test distance_to_in(shell, Point3{T}(0, 0, 0), vx) ≈ 40         # from the cavity
        @test distance_to_in(shell, Point3{T}(-100, 0, 0), vx) ≈ 50
        @test distance_to_in(upper, Point3{T}(0, 0, -100), vz) ≈ 140     # through the cavity
        @test distance_to_in(upper, Point3{T}(45, 0, -10), vz) ≈ 10      # through the z = 0 face
        @test distance_to_in(seg, Point3{T}(-45, 45, 0), vx) ≈ 45        # through the ϕ = π/2 face
        @test distance_to_out(orb, Point3{T}(0, 0, 0), vx) ≈ 50
        @test distance_to_out(shell, Point3{T}(-45, 0, 0), vx) ≈ 5       # into the cavity
        @test distance_to_out(upper, Point3{T}(0, 45, 1), -vz) ≈ 1       # through the z = 0 face
        @test distance_to_out(seg, Point3{T}(1, 45, 0), -vx) ≈ 1         # through the ϕ = π/2 face
        @test distance_to_in(cap, Point3{T}(-100, 0, 30), vx) ≈ 100 - 30tan(π/6)
        @test distance_to_out(cap, Point3{T}(0, 0, 30), vx) ≈ 30tan(π/6)
        @test distance_to_in(cap, Point3{T}(0, 0, -100), vz) ≈ 100       # through the apex
        @test distance_to_in(band, Point3{T}(30, 0, -100), vz) ≈ 70       # through the 135° cone
        @test distance_to_out(band, Point3{T}(30, 0, 0), vz) ≈ 30         # through the 45° cone
        @test distance_to_in(band, Point3{T}(0, 0, 0), vx) ≈ 20           # through the inner sphere
        @test distance_to_in(band, Point3{T}(5, 0, 100), -vz) == Inf      # along the axis, in the cones
        @test distance_to_out(wide, Point3{T}(10, 0, -1), -vz) ≈ 10tan(π/6) - 1 rtol=1e-3
        # entry through a ϕ plane, at x = 0: the other ϕ plane is behind (concave wedge)
        @test distance_to_in(upcone, Point3{T}(38, -16, 42), normalize(Vector3{T}(-38, -4, -7.5))) ≈ 38 * sqrt(1 + (4/38)^2 + (7.5/38)^2) rtol=1e-4
        @test distance_to_out(orb, Point3{T}(60, 0, 0), vx) == -1
        @test distance_to_in(orb, Point3{T}(50, 0, 0), -vx) <= 0         # on the surface, moving in
        @test distance_to_in(orb, Point3{T}(50, 0, 0), vx) == Inf        # on the surface, moving out
        @test distance_to_out(orb, Point3{T}(50, 0, 0), vx) ≈ 0 atol=tol
    end

    @testset "normals" begin
        @test normal(orb, Point3{T}(50, 0, 0)) ≈ vx
        @test normal(shell, Point3{T}(40, 0, 0)) ≈ -vx
        @test normal(upper, Point3{T}(45, 0, 0)) ≈ -vz
        @test normal(seg, Point3{T}(45, 0, 10)) ≈ -vy
        @test normal(seg, Point3{T}(0, 45, 10)) ≈ -vx
        @test normal(orb, Point3{T}(20, 0, 0)) == zero(Vector3{T})
        @test normal(cap, Point3{T}(20, 0, 40cos(π/6))) ≈ Vector3{T}(cos(π/6), 0, -sin(π/6))
        @test normal(band, Point3{T}(0, 30, -30)) ≈ Vector3{T}(0, -sqrt(T(0.5)), -sqrt(T(0.5)))
    end

    @testset "capacity and surface" begin
        @test capacity(orb) ≈ 4π/3 * 50^3
        @test capacity(upper) ≈ 2π/3 * (50^3 - 40^3)
        @test surface(orb) ≈ 4π * 50^2
        @test capacity(band) ≈ 2π/3 * (50^3 - 20^3) * 2cos(π/4)
        @test surface(cap) ≈ 50^2 * 2π * (1 - cos(π/6)) + 50^2 * sin(π/6) * π
        for s in (shell, upper, seg, seg2, cap, band, low, wide, upcone)
            @test estimate_capacity(s; npoints=200_000) ≈ capacity(s) rtol=0.03
        end
    end

    # Random rays: the distances must lead to points on the surface, and a ray that enters the shape
    # must then be inside it
    @testset "random rays $(s)" for s in (orb, shell, upper, lower, seg, seg2, cap, band, low, wide, upcone)
        rng = MersenneTwister(1234)
        nin = nout = okin = okout = nhit = 0
        for _ in 1:4000
            p = Point3{T}((rand(rng, T, 3) .* 2 .- 1) .* 80)
            d = normalize(Vector3{T}(randn(rng, T, 3)))
            pos = inside(s, p)
            if pos == kOutside
                nin += 1
                t = distance_to_in(s, p, d)
                if t == Inf
                    okin += 1
                else
                    nhit += 1
                    q = p + t * d
                    okin += inside(s, q) == kSurface && inside(s, q + 10tol * d) != kOutside
                end
            elseif pos == kInside
                nout += 1
                t = distance_to_out(s, p, d)
                q = p + t * d
                okout += inside(s, q) == kSurface && inside(s, q + 10tol * d) != kInside
            end
        end
        @test nhit > 100
        @test okin >= 0.995 * nin            # a few rays may graze the surfaces
        @test okout >= 0.995 * nout
    end
end
