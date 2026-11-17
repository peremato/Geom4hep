using Random

# tree volumes in the order used by close_geometry (preorder, shared volumes once)
function volume_order(world)
    vols = Any[]; seen = IdDict{Any,Bool}()
    walk(v) = (haskey(seen, v) && return; seen[v] = true; push!(vols, v); foreach(d -> walk(d.volume), v.daughters))
    walk(world)
    return vols
end

@testset "GeometryModel{$T}" for T in (Float64, Float32)
    P(x, y, z) = Point3{T}(x, y, z)
    tr(x, y, z) = Transformation3D{T}(x, y, z)
    mat = Material{T}("mat"; density=1)

    @testset "boolean shapes as CSG" begin
        b = Box{T}(1, 1, 1)
        chain = BooleanUnion(BooleanUnion(BooleanUnion(b, b, tr(2, 0, 0)), b, tr(4, 0, 0)), b, tr(6, 0, 0))
        sub = BooleanSubtraction(BooleanSubtraction(Box{T}(10, 10, 10), b, tr(2, 0, 0)), b, tr(-2, 0, 0))
        world = Volume{T}("world", Box{T}(100, 100, 100), mat)
        place_daughter!(world, tr(0, 0, 0), Volume{T}("chain", chain, mat))
        place_daughter!(world, tr(0, 50, 0), Volume{T}("sub", sub, mat))
        m = close_geometry(world)
        kinds = eltype.(m.shapes)
        @test Geom4hep.FlatCSG{T} in kinds
        @test count(==(Geom4hep.FlatCSG{T}), kinds) == 1
        @test length(m.shapes[findfirst(==(Geom4hep.FlatCSG{T}), kinds)]) == 2     # two boolean shapes
        @test length(m.csg_components) == 4 + 3                    # 4 union components + box minus 2 boxes
        @test length(m.csg_program) == (4 + 3) + (3 + 2)            # components + operations
        @test Box{T} in kinds && length(m.shapes[findfirst(==(Box{T}), kinds)]) == 3   # world, b, big box
        rng = MersenneTwister(1)
        for (k, v) in enumerate(volume_order(world))
            lo, hi = extent(v.shape)
            for _ in 1:500
                p = Point3{T}(lo .- 1 .+ (hi .- lo .+ 2) .* rand(rng, T, 3))
                @test inside(v.shape, p) == shape_call(inside, m, m.volumes[k].shape, p)
            end
        end
    end

    @testset "two-component booleans" begin
        # flat booleans of two components (fast path of distance_to_out for A - B and A ∩ B) against the tree shapes
        rot(x, y, z, a, b, c) = Transformation3D{T}(x, y, z, RotXYZ{T}(a, b, c))
        shapes = [BooleanSubtraction(Box{T}(350, 350, 25), Tube{T}(0, 32, 50, 0, 2π), tr(0, 0, 0)),
                  BooleanSubtraction(Tube{T}(5, 20, 20, 0, 2π), Box{T}(8, 8, 30), rot(10, 0, 0, 0, 0, 0.5)),
                  BooleanSubtraction(Sphere{T}(0, 20, 0, 2π, 0, π), Box{T}(10, 10, 10), tr(15, 0, 0)),
                  BooleanIntersection(Box{T}(20, 20, 20), Sphere{T}(10, 25, 0, 2π, 0, π), tr(0, 0, 5)),
                  BooleanIntersection(Tube{T}(0, 20, 20, 0, 2π), Box{T}(15, 15, 15), rot(3, 0, 0, 0.4, 0, 0.7)),
                  BooleanUnion(Box{T}(10, 10, 10), Tube{T}(0, 5, 20, 0, 2π), tr(0, 0, 0))]
        world = Volume{T}("world", Box{T}(1000, 1000, 1000), mat)
        foreach(s -> place_daughter!(world, tr(0, 0, 0), Volume{T}("s", s, mat)), shapes)
        m = close_geometry(world)
        rng = MersenneTwister(3)
        nbad = ninside = nnormal = 0
        for (k, s) in enumerate(shapes), _ in 1:2000
            ref = m.volumes[k+1].shape
            lo, hi = extent(s)
            p = Point3{T}(lo .- 2 .+ (hi .- lo .+ 4) .* rand(rng, T, 3))
            d = Vector3{T}(normalize(randn(rng, 3))...)
            pos = inside(s, p)
            nbad += pos != shape_call(inside, m, ref, p)
            pos == kInside || continue
            ninside += 1
            dist = distance_to_out(s, p, d)
            nbad += !isapprox(dist, shape_call(distance_to_out, m, ref, p, d); atol=10push_tolerance(T), rtol=1e-4)
            n = shape_call(normal, m, ref, p + dist * d)          # at the exit point, along the ray
            nnormal += isapprox(norm(n), 1; rtol=1e-3) && n ⋅ d > 0
        end
        @test nbad == 0
        @test nnormal == ninside
    end

    @testset "GDML models" begin
        for f in ("boxes.gdml", "polycones.gdml", "cuttubes.gdml", "booleans.gdml")
            world = read_gdml(joinpath(@__DIR__, "..", "examples", f), T)
            m = close_geometry(world)
            vols = volume_order(world)
            @test length(m.volumes) == length(vols)
            @test length(m.placements) == sum(v -> length(v.daughters), vols)
            rng = MersenneTwister(2)
            nbad = 0
            for (k, v) in enumerate(vols), _ in 1:100
                lo, hi = extent(v.shape)
                p = Point3{T}(lo .+ (hi .- lo) .* rand(rng, T, 3))
                d = Vector3{T}(normalize(randn(rng, 3))...)
                ref = m.volumes[k].shape
                pos = inside(v.shape, p)
                nbad += pos != shape_call(inside, m, ref, p)
                if pos == kInside
                    nbad += !isapprox(distance_to_out(v.shape, p, d), shape_call(distance_to_out, m, ref, p, d); atol=10push_tolerance(T), rtol=1e-4)
                elseif pos == kOutside
                    nbad += !isapprox(distance_to_in(v.shape, p, d), shape_call(distance_to_in, m, ref, p, d); atol=10push_tolerance(T), rtol=1e-4)
                end
            end
            @test nbad == 0
        end
    end

    @testset "polycone views" begin
        world = read_gdml(joinpath(@__DIR__, "..", "examples", "polycones.gdml"), T)
        m = close_geometry(world)
        @test !isempty(m.cones)
        @test Geom4hep.FlatPolycone{T} in eltype.(m.shapes)
        @test all(isbitstype, eltype.(m.shapes))
    end
end

@testset "Navigator{$T}" for T in (Float64, Float32)
    mat   = Material{T}("vacuum"; density=0)
    world = Volume{T}("world", Box{T}(100, 100, 100), mat)
    inner = Volume{T}("inner", Box{T}(10, 10, 10), Material{T}("iron"; density=7.87))
    core  = Volume{T}("core", Box{T}(2, 2, 2), Material{T}("lead"; density=11.35))
    place_daughter!(inner, Transformation3D{T}(0, 0, 0), core)
    place_daughter!(world, Transformation3D{T}(0, 0, 0), inner)
    place_daughter!(world, Transformation3D{T}(50, 0, 0), inner)
    m = close_geometry(world; bvh_threshold=1)          # also exercise the BVH with few daughters

    @test locate(m, Point3{T}(0, 0, 0)).volume == 3     # world=1, inner=2, core=3 (preorder)
    @test locate(m, Point3{T}(55, 0, 0)).volume == 2
    @test locate(m, Point3{T}(30, 0, 0)).volume == 1
    @test !in_world(locate(m, Point3{T}(300, 0, 0)))
    @test locate(m, Point3{T}(51, 1, 1)).depth == 2

    point, dir = Point3{T}(-90, 0, 0), Vector3{T}(1, 0, 0)
    s = locate(m, point)
    steps, vols = T[], Int32[]
    while in_world(s) && length(steps) < 50
        step, s = compute_step(m, s, point, dir, T(1e6))
        push!(steps, step)
        point += step * dir
        in_world(s) && push!(vols, s.volume)
    end
    @test steps ≈ T[80, 8, 4, 8, 30, 8, 4, 8, 40] atol = 10push_tolerance(T)
    @test vols == Int32[2, 3, 2, 1, 2, 3, 2, 1]

    # step limited by step_limit: same state
    s0 = locate(m, Point3{T}(-90, 0, 0))
    step, s1 = compute_step(m, s0, Point3{T}(-90, 0, 0), dir, T(5))
    @test step == 5 && s1 == s0

    # X-ray: rays along z through the world; pixels on the axis x=y=0 cross 20 iron + 4 lead (- 4 iron)
    u, v, img, ns = xray(m, 3, 201)
    i, j = argmin(abs.(u)), argmin(abs.(v))
    @test img[i, j] ≈ 16 * T(7.87) + 4 * T(11.35) rtol = 1e-3
    @test maximum(img) ≈ img[i, j] rtol = 1e-3
    @test minimum(img) == 0
end

@testset "Navigation at surfaces{$T}" for T in (Float64, Float32)
    # a daughter sharing the inner surface of its mother, crossed by rays at grazing angles: the
    # point pushed into the mother is still on the surface of the daughter, which must be entered
    world  = Volume{T}("world", Box{T}(5000, 5000, 100), Material{T}("vacuum"; density=0))
    mother = Volume{T}("mother", Tube{T}(3800, 4000, 50, 0, 2π), Material{T}("air"; density=1))
    inner  = Volume{T}("inner", Tube{T}(3800, 3900, 50, 0, 2π), Material{T}("iron"; density=10))
    place_daughter!(mother, Transformation3D{T}(0, 0, 0), inner)
    place_daughter!(world, Transformation3D{T}(0, 0, 0), mother)
    m = close_geometry(world)
    nbad = 0
    for y in range(T(3700), T(3799), length=200)
        p, d = Point3{T}(-4500, y, 0), Vector3{T}(1, 0, 0)
        s = locate(m, p)
        mass = zero(T)
        while in_world(s)
            ρ = m.materials[material_index(m, s)].density
            step, s = compute_step(m, s, p, d, T(Inf))
            p += step * d
            mass += step * ρ
        end
        x(r) = sqrt(r^2 - y^2)
        expected = 2 * (10 * (x(T(3900)) - x(T(3800))) + (x(T(4000)) - x(T(3900))))
        nbad += !isapprox(mass, expected; rtol=1e-3)
    end
    @test nbad == 0

    # the same with boxes (whose distance_to_in on the surface is not exactly 0): a daughter sharing
    # the +y face of its mother, rays entering the mother through that face at grazing angles
    world  = Volume{T}("world", Box{T}(500, 500, 500), Material{T}("vacuum"; density=0))
    mother = Volume{T}("mother", Box{T}(100, 100, 100), Material{T}("air"; density=1))
    inner  = Volume{T}("inner", Box{T}(100, 10, 100), Material{T}("iron"; density=10))
    place_daughter!(mother, Transformation3D{T}(0, 90, 0), inner)
    place_daughter!(world, Transformation3D{T}(0, 0, 0), mother)
    m = close_geometry(world)
    nbad = 0
    for θ in range(T(0.01), T(0.2), length=50), z in (T(-50), T(0), T(50))
        d = Vector3{T}(cos(θ), -sin(θ), 0)
        p = Point3{T}(-50, 100, z) - T(200) * d                # crossing y = 100 at x = -50
        s = locate(m, p)
        mass = zero(T)
        while in_world(s)
            ρ = m.materials[material_index(m, s)].density
            step, s = compute_step(m, s, p, d, T(Inf))
            p += step * d
            mass += step * ρ
        end
        l, L = 20 / sin(θ), 150 / cos(θ)                       # paths to the bottom of the daughter, to x = 100
        expected = l < L ? 10l + (L - l) : 10L
        nbad += !isapprox(mass, expected; rtol=1e-3)
    end
    @test nbad == 0

    # leaving a (non-convex) boolean shape at an angle from its surface: not entering it again there
    u = BooleanUnion(Box{T}(10, 10, 10), Box{T}(10, 50, 50), Transformation3D{T}(30, 0, 0))
    world = Volume{T}("world", Box{T}(100, 100, 100), Material{T}("vacuum"; density=0))
    place_daughter!(world, Transformation3D{T}(0, 0, 0), Volume{T}("u", u, Material{T}("mat"; density=1)))
    m = close_geometry(world)
    nbad = 0
    for y in range(T(-9), T(9), length=50), z in range(T(-9), T(9), length=50)
        p = Point3{T}(10, y, z)                            # on the +x face of the first box
        d = Vector3{T}(1 / 2, 0, sqrt(T(3)) / 2)           # 60° from the normal
        nbad += !isapprox(shape_call(distance_to_in, m, m.volumes[2].shape, p, d), 20; atol=1e-2)
    end
    @test nbad == 0
end

@testset "Navigator agreement (GDML, flat model)" begin
    T = Float64
    world = read_gdml(joinpath(@__DIR__, "..", "examples", "boxes.gdml"), T)
    m = close_geometry(world)
    old = NavigatorState{T}(world, TrivialNavigator(world))
    lo, hi = extent(world.shape)
    rng = MersenneTwister(1234)
    for _ in 1:200
        p0 = Point3{T}(lo .+ (hi .- lo) .* (T(0.1) .+ T(0.8) .* rand(rng, T, 3)))
        d  = Vector3{T}(normalize(randn(rng, T, 3)))
        s = locate(m, p0); locate_global_point!(old, p0)
        pa, pb = p0, p0
        while in_world(s)
            a, s = compute_step(m, s, pa, d, T(1e10))
            b = compute_step!(old, pb, d, T(1e10))
            @test a ≈ b
            pa += a * d; pb += b * d
        end
        @test !is_in_volume(old)
    end
end

@testset "Navigation robustness (cms2018)" begin
    # rays from the origin used to fail with the tree navigator ("Negative distanceToOut")
    T = Float64
    world = read_gdml(joinpath(@__DIR__, "..", "examples", "cms2018.gdml"), T)
    m = close_geometry(world)
    rng = MersenneTwister(1)
    for _ in 1:20
        d = Vector3{T}(normalize(randn(rng, 3))...)
        p = Point3{T}(0, 0, 0)
        s = locate(m, p)
        n = 0
        while in_world(s) && n < 20_000
            step, s = compute_step(m, s, p, d, T(1e10))
            p += step * d
            n += 1
        end
        @test !in_world(s)
        @test n < 20_000
    end
end
