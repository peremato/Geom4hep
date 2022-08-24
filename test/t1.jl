#---Testing-----------------------------------------------------------------------------------------
using Revise
using Geom4hep
using GeometryBasics

const T = Float64

box = Box{T}(1,1,1)
cone = Cone{T}(0,0.5, 0, 0, 0.5, 0, 2π)
tube = Tube{T}(0,0.5,0.5, 0,2π)
mbox = Box{T}(0.5, 0.5, 0.5)
water = Material{T}("water"; density=1.)
vacuum = Material{T}("vacuum"; density=0 )
#vbox = Volume{T}("aggregate", Box{T}(2,2,2), vacuum)
vbox = Assembly{T}("aggregate")
placeDaughter!(vbox, one(Transformation3D{T}), Volume{T}("box", box, water))
placeDaughter!(vbox, Transformation3D{T}(1.5,0,0, 0, -π/2, 0), Volume{T}("cone",cone, water))
placeDaughter!(vbox, Transformation3D{T}(0,1.5,0, π/2, 0, 0), Volume{T}("tube",tube, water))
placeDaughter!(vbox, Transformation3D{T}(0,0,1.5), Volume{T}("mbox",mbox, water))

assembly = Assembly{T}("myassembly")
placeDaughter!(assembly, Transformation3D{T}(2,0,2), vbox)
placeDaughter!(assembly, Transformation3D{T}(2,0,-2, π/4, 0, 0), vbox)
placeDaughter!(assembly, Transformation3D{T}(-2,0,-2, 0, π/4, 0), vbox)
placeDaughter!(assembly, Transformation3D{T}(-2,0,2, 0, 0, π/4), vbox)

world = Volume{T}("world", Box{T}(4,4,4), vacuum)
placeDaughter!(world, Transformation3D{T}(0, 2,0), assembly)
placeDaughter!(world, Transformation3D{T}(0,-2,0), assembly)
#placeDaughter!(world, Transformation3D{T}(0,0,0, 0,0,π/4), vbox)
draw(world)

using GLMakie

rx = generateXRay(world, 1e6, 1)
ry = generateXRay(world, 1e6, 2)
rz = generateXRay(world, 1e6, 3)
#----Plot the results--------------------------------------
fig = Figure(resolution = (1000, 1000))
heatmap!(Axis(fig[1, 1], title = "X direction"), rx..., colormap=:grayC, colorrange=(0,maximum(rx[3])))
heatmap!(Axis(fig[2, 1], title = "Y direction"), ry..., colormap=:grayC, colorrange=(0,maximum(ry[3])))
heatmap!(Axis(fig[1, 2], title = "Z direction"), rz..., colormap=:grayC, colorrange=(0,maximum(rz[3])))
draw(LScene(fig[2, 2]), world, 3)
display(fig)


#----Testing BVH---------------------------------------------
using Revise
using Geom4hep
using GeometryBasics

const T = Float64

world = processGDML("examples/trackML.gdml", Float64)
vol = world.daughters[2].volume.daughters[1].volume.daughters[1].volume

using GLMakie
draw(vol, wireframe=true, bvh=true)

nav = BVHNavigator{T}(world)
state1 = NavigatorState{T}(world)
state2 = NavigatorState{T}(world, nav)

p = Point3{T}(25,25,0)
using BenchmarkTools
@btime locateGlobalPoint!(state1, p)
@btime locateGlobalPoint!(state2, p)



