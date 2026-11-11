using Random

@testset "BVH{$T}" for T in (Float64, Float32)
    mat = Material{T}("vacuum"; density=0)
    world = Volume{T}("world", Box{T}(1000, 1000, 1000), mat)
    cell = Volume{T}("cell", Box{T}(4, 4, 4), mat)
    for i in -5:5, j in -5:5
        placeDaughter!(world, Transformation3D{T}(10i, 10j, 0), cell)
    end
    n = length(world.daughters)
    bvh = buildBVH(world.daughters)
    @test sort(pvolindices(bvh)) == 1:n
    # Only the daughter at the origin is selected when the AABB filter contains the origin
    p = Point3{T}(1, 1, 1)
    sel = pvolindices(x -> Geom4hep.inside(x, p), bvh)
    @test any(world.daughters[i].transformation.translation ≈ Vector3{T}(0, 0, 0) for i in sel)
    @test length(sel) < n
end

@testset "Navigation{$T}" for T in (Float64, Float32)
    mat   = Material{T}("vacuum"; density=0)
    world = Volume{T}("world", Box{T}(100, 100, 100), mat)
    inner = Volume{T}("inner", Box{T}(10, 10, 10), mat)
    core  = Volume{T}("core", Box{T}(2, 2, 2), mat)
    placeDaughter!(inner, Transformation3D{T}(0, 0, 0), core)
    placeDaughter!(world, Transformation3D{T}(0, 0, 0), inner)
    placeDaughter!(world, Transformation3D{T}(50, 0, 0), inner)   # inner centred at x = +50

    for nav in (TrivialNavigator(world), BVHNavigator(world))
        state = NavigatorState{T}(world, nav)
        @test locateGlobalPoint!(state, Point3{T}(0, 0, 0))
        @test currentVolume(state) === core
        @test locateGlobalPoint!(state, Point3{T}(55, 0, 0))
        @test currentVolume(state) === inner
        @test locateGlobalPoint!(state, Point3{T}(30, 0, 0))
        @test currentVolume(state) === world
        @test !locateGlobalPoint!(state, Point3{T}(300, 0, 0))

        # Shoot along +x from x = -90 and record the boundaries crossed
        point = Point3{T}(-90, 0, 0)
        dir   = Vector3{T}(1, 0, 0)
        locateGlobalPoint!(state, point)
        steps, labels = T[], String[]
        while isInVolume(state) && length(steps) < 50
            step = computeStep!(state, point, dir, T(1e6))
            push!(steps, step)
            point = point + dir * step
            isInVolume(state) && push!(labels, currentVolume(state).label)
        end
        expected = T[80, 8, 4, 8, 30, 8, 4, 8, 40]   # boundaries at -10,-2,2,10,40,48,52,60,100
        @test length(steps) == length(expected)
        @test all(isapprox.(steps, expected; atol=10kPushTolerance(T)))
        @test labels == ["inner", "core", "inner", "world", "inner", "core", "inner", "world"]
    end
end

@testset "Navigator agreement (GDML)" begin
    T = Float64
    world = processGDML(joinpath(@__DIR__, "..", "examples", "boxes.gdml"), T)
    triv  = NavigatorState{T}(world, TrivialNavigator(world))
    bvh   = NavigatorState{T}(world, BVHNavigator(world))
    lo, hi = extent(world.shape)
    rng = MersenneTwister(1234)
    ntotal = 0
    for _ in 1:200
        p0 = Point3{T}(lo .+ (hi .- lo) .* (T(0.1) .+ T(0.8) .* rand(rng, T, 3)))
        d  = Vector3{T}(normalize(randn(rng, T, 3)))
        pa, pb = p0, p0
        locateGlobalPoint!(triv, pa); locateGlobalPoint!(bvh, pb)
        @test currentVolume(triv) === currentVolume(bvh)
        n = 0
        while isInVolume(triv) && n < 1000
            sa = computeStep!(triv, pa, d, T(1e10))
            sb = computeStep!(bvh,  pb, d, T(1e10))
            sa ≈ sb || break
            pa += d * sa; pb += d * sb; n += 1
        end
        ntotal += n
        @test isInVolume(triv) == isInVolume(bvh)
        @test pa ≈ pb
    end
    @test ntotal > 200
end
