@testset "GDML{$T}" for T in (Float64, Float32)
    world = read_gdml(joinpath(@__DIR__, "twoboxes.gdml"), T)
    @test world isa Volume{T}
    @test world.label == "World"
    @test world.shape isa Box{T}
    @test extent(world.shape)[2] ≈ Point3{T}(500, 500, 500)     # 100 cm full length → 500 mm half
    @test world.material.name == "Air"
    @test world.material.density ≈ T(0.00121399936124299)
    @test length(world.material.composition) == 4
    @test world.material.composition[1].element.name == "Nitrogen"
    @test length(world.material.composition[1].element.composition) == 2

    @test length(world.daughters) == 1
    trd = world.daughters[1]
    @test trd.volume.label == "TRD"
    @test trd.volume.shape isa Trd{T}
    @test trd.volume.material.name == "Steel"
    # position (0,0,200) mm, rotation 45° around z
    @test trd.transformation.translation ≈ Vector3{T}(0, 0, 200)
    @test Geom4hep.contains(trd, Point3{T}(0, 0, 200))
    @test !Geom4hep.contains(trd, Point3{T}(0, 0, 0))

    @testset "sphere and optical elements" begin
        w = @test_logs read_gdml(joinpath(@__DIR__, "mirror.gdml"), T)      # no warnings
        mirror = w.daughters[1].volume
        @test mirror.shape isa BooleanIntersection
        @test mirror.shape.left isa Sphere{T}
        model, index = close_geometry(w; with_index=true)
        @test [v.label for v in index.volumes] == ["World", "Mirror"]
        @test [m.name for m in index.materials] == ["Gas", "Glass"]
        # a ray along z from the origin hits the inner surface of the mirror at z = 1500 - 500
        p, d = Point3{T}(0, 0, 0), Vector3{T}(0, 0, 1)
        s = locate(model, p)
        @test index.volumes[s.volume].label == "World"
        step, s = compute_step(model, s, p, d, T(Inf))
        @test step ≈ 1000 atol = 10Geom4hep.tolerance(T)
        @test index.volumes[s.volume].label == "Mirror"
        @test index.materials[material_index(model, s)].name == "Glass"
    end

    @testset "$f" for f in ("boxes.gdml", "booleans.gdml", "cuttubes.gdml", "polycones.gdml")
        w = read_gdml(joinpath(@__DIR__, "..", "examples", f), T)
        @test w isa Volume{T}
        if f == "booleans.gdml"     # a single boolean solid as world
            @test w.shape isa Geom4hep.BooleanShape
        else
            @test !isempty(w.daughters)
        end
    end
end

@testset "GDML units" begin
    @test Geom4hep._lunit(Dict("lunit" => "cm")) == 10
    @test Geom4hep._lunit(Dict{String,String}()) == 1               # default mm
    @test Geom4hep._aunit(Dict("aunit" => "deg")) ≈ π/180
    @test Geom4hep._aunit(Dict{String,String}()) == 1               # default rad
    @test_throws ArgumentError Geom4hep._lunit(Dict("lunit" => "furlong"))
end
