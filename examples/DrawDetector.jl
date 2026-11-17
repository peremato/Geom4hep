
using GeometryBasics
using Geom4hep: Volume, Box, Sphere, Trd, Material, Transformation3D, RotXYZ, place_daughter!, draw!

# Geometry construction
function build_geom(T::Type) 
    world = Volume{T}("World", Box{T}(100,100,100), Material{T}("vacuum"; density=0.0))
    box1 = Volume{T}("box1", Box{T}(10,20,30), Material{T}("iron"; density=7.0))
    box2 = Volume{T}("box2", Box{T}(4,4,4), Material{T}("gold"; density=19.0))
    place_daughter!(box1, Transformation3D{T}(0,0,0),box2)
    place_daughter!(world, Transformation3D{T}( 50, 50, 50, RotXYZ{T}(π/4, 0, 0)), box1)
    place_daughter!(world, Transformation3D{T}( 50, 50,-50, RotXYZ{T}(0, π/4, 0)), box1)
    place_daughter!(world, Transformation3D{T}( 50,-50, 50, RotXYZ{T}(0, 0, π/4)), box1)
    place_daughter!(world, Transformation3D{T}( 50,-50,-50, RotXYZ{T}(π/4, 0, 0)), box1)
    place_daughter!(world, Transformation3D{T}(-50, 50, 50, RotXYZ{T}(0, π/4, 0)), box1)
    place_daughter!(world, Transformation3D{T}(-50,-50, 50, RotXYZ{T}(0, 0, π/4)), box1)
    place_daughter!(world, Transformation3D{T}(-50, 50,-50, RotXYZ{T}(π/4, 0, 0)), box1)
    place_daughter!(world, Transformation3D{T}(-50,-50,-50, RotXYZ{T}(0, π/4, 0)), box1)
    trd1 = Volume{T}("trd1", Trd{T}(25,10,25,10,10), Material{T}("water", density=1.0))
    place_daughter!(world, Transformation3D{T}(0,0,0), trd1)
    ball = Volume{T}("ball", Sphere{T}(0, 10), Material{T}("lead", density=11.0))
    place_daughter!(trd1, Transformation3D{T}(0,0,0), ball)
    return world
end

world = build_geom(Float64)
# or 
#world = read_gdml("examples/boxes.gdml")

using GLMakie
fig = Figure()
scene = LScene(fig[1, 1])
draw!(scene, world)
display(fig)



