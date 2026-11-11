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
