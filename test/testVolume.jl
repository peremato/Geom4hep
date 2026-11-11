@testset "Volume{$T}" for T in (Float64, Float32)
    air   = Material{T}("Air"; density=0.0012)
    lead  = Material{T}("Lead"; density=11.35)
    world = Volume{T}("world", Box{T}(100, 100, 100), air)
    small = Volume{T}("small", Box{T}(5, 5, 5), lead)
    placeDaughter!(world, Transformation3D{T}(50, 0, 0), small)
    placeDaughter!(world, Transformation3D{T}(-50, 0, 0), small)

    @test length(world.daughters) == 2
    @test world[1] === small
    @test world.id != small.id
    @test children(world) == [small, small]
    @test Geom4hep.contains(world.daughters[1], Point3{T}(52, 0, 0))
    @test !Geom4hep.contains(world.daughters[1], Point3{T}(-52, 0, 0))
    bb = world.daughters[2].aabb
    @test bb.min ≈ Point3{T}(-55, -5, -5)
    @test bb.max ≈ Point3{T}(-45, 5, 5)

    @testset "Assembly" begin
        asm = Assembly{T}("asm")
        placeDaughter!(asm, Transformation3D{T}(0, 20, 0), small)
        placeDaughter!(asm, Transformation3D{T}(0, -20, 30), small)
        lo, hi = extent(asm.shape)
        @test lo ≈ Point3{T}(-5, -25, -5)
        @test hi ≈ Point3{T}(5, 25, 35)
        @test inside(asm.shape, Point3{T}(0, 20, 0)) == kInside
        @test inside(asm.shape, Point3{T}(0, 0, 0)) == kOutside
        # placing an assembly flattens its daughters into the mother
        top = Volume{T}("top", Box{T}(100, 100, 100), air)
        placeDaughter!(top, Transformation3D{T}(10, 0, 0), asm)
        @test length(top.daughters) == 2
        @test all(d.volume === small for d in top.daughters)
    end

    @testset "getWorld" begin
        w = getWorld(small)
        @test w.label == "world"
        @test length(w.daughters) == 1
        @test getWorld(world) === world
    end
end

@testset "Tolerances" begin
    @test kTolerance() == kTolerance(Float64)
    @test kPushTolerance() == kPushTolerance(Float64)
    @test kPushTolerance(Float32) isa Float32
end

@testset "Mesh" begin
    T = Float64
    shapes = [Box{T}(1, 2, 3), Tube{T}(1, 2, 3, 0, 2π), Tube{T}(1, 2, 3, 0, π), Trd{T}(1, 2, 3, 4, 5),
              Geom4hep.Cone{T}(1, 2, 1, 3, 4, 0, 2π), Trap{T}(2, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0),
              Polycone{T}([0., 1, 1], [2., 3, 2], [-1., 0, 1], 0, 2π),
              CutTube{T}(1, 2, 3, 0, 2π, Vector3{T}(0, -0.7, -0.71), Vector3{T}(0.7, 0, 0.71)),
              BooleanUnion(Box{T}(1, 1, 1), Box{T}(1, 1, 1), Transformation3D{T}(1, 0, 0))]
    for s in shapes
        m = mesh(s)
        @test length(coordinates(m)) > 0
        @test length(faces(m)) > 0
    end
end
