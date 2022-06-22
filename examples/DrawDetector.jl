
using Revise
using GeometryBasics
using Geom4hep
using GLMakie

# Geometry construction
function buildGeom(T::Type) 
    world = Volume{T}("World", Box{T}(100,100,100), Material{T}("vacuum"; density=0.0))
    box1 = Volume{T}("box1", Box{T}(10,20,30), Material{T}("iron"; density=7.0))
    box2 = Volume{T}("box2", Box{T}(4,4,4), Material{T}("gold"; density=19.0))
    placeDaughter!(box1, Transformation3D{T}(0,0,0),box2)
    placeDaughter!(world, Transformation3D{T}( 50, 50, 50, RotXYZ{T}(π/4, 0, 0)), box1)
    placeDaughter!(world, Transformation3D{T}( 50, 50,-50, RotXYZ{T}(0, π/4, 0)), box1)
    placeDaughter!(world, Transformation3D{T}( 50,-50, 50, RotXYZ{T}(0, 0, π/4)), box1)
    placeDaughter!(world, Transformation3D{T}( 50,-50,-50, RotXYZ{T}(π/4, 0, 0)), box1)
    placeDaughter!(world, Transformation3D{T}(-50, 50, 50, RotXYZ{T}(0, π/4, 0)), box1)
    placeDaughter!(world, Transformation3D{T}(-50,-50, 50, RotXYZ{T}(0, 0, π/4)), box1)
    placeDaughter!(world, Transformation3D{T}(-50, 50,-50, RotXYZ{T}(π/4, 0, 0)), box1)
    placeDaughter!(world, Transformation3D{T}(-50,-50,-50, RotXYZ{T}(0, π/4, 0)), box1)
    trd1 = Volume{T}("trd1", Trd{T}(25,10,25,10,20), Material{T}("water", density=1.0))
    placeDaughter!(world, Transformation3D{T}(0,0,0), trd1)
    return world
end

world = buildGeom(Float64)
# or 
world = processGDML("examples/boxes.gdml")

fig = Figure()
scene = LScene(fig[1, 1])
draw(scene, world)
display(fig)



