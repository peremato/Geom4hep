using Geom4hep, Geom4hep.Units
using CairoMakie: CairoMakie, Figure, LScene

const T = Float64;

vacuum = Material{T}("vacuum"; density=0)
air    = Material{T}("air"; density=0.0012)
iron   = Material{T}("iron"; density=7.87)
lead   = Material{T}("lead"; density=11.35)
water  = Material{T}("water"; density=1.0);

world  = Volume{T}("world", Box{T}(1m, 1m, 1m), air)
plate  = Volume{T}("plate", Box{T}(10cm, 20cm, 30cm), iron)
brick  = Volume{T}("brick", Box{T}(4cm, 4cm, 4cm), lead)
tank   = Volume{T}("tank", Trd{T}(25cm, 10cm, 25cm, 10cm, 10cm), water)
ball   = Volume{T}("ball", Sphere{T}(0, 8cm), lead)
pipe   = Volume{T}("pipe", Tube{T}(4cm, 5cm, 80cm, 0, 2π), iron);

place_daughter!(plate, Transformation3D{T}(0, 0, 0), brick)
for (i, (x, y, z)) in enumerate(Iterators.product((-50cm, 50cm), (-50cm, 50cm), (-50cm, 50cm)))
    rotation = RotXYZ{T}(i % 3 == 0 ? π/4 : 0, i % 3 == 1 ? π/4 : 0, i % 3 == 2 ? π/4 : 0)
    place_daughter!(world, Transformation3D{T}(x, y, z, rotation), plate)
end
place_daughter!(tank, Transformation3D{T}(0, 0, 0), ball)
place_daughter!(world, Transformation3D{T}(0, 0, 0), tank)
place_daughter!(world, Transformation3D{T}(0, 0, 40cm, π/2, 0, 0), pipe) # along y, above the tank

world

world[1, 1]

bar = Volume{T}("bar", Box{T}(20cm, 2cm, 2cm), iron)
group = Assembly{T}("bars")
place_daughter!(group, Transformation3D{T}(0, -5cm, 0), bar)
place_daughter!(group, Transformation3D{T}(0, 5cm, 0), bar)
place_daughter!(world, Transformation3D{T}(0, 0, -30cm), group)
length(world.daughters)

fig = Figure(size=(800, 600))
draw!(LScene(fig[1, 1], show_axis=false), world)
fig

fig = Figure(size=(800, 600))
draw!(LScene(fig[1, 1], show_axis=false), world; wireframe=true, maxlevel=2)
fig
