using Geom4hep
using CairoMakie: CairoMakie, Figure, LScene, Axis, heatmap!, DataAspect

gdmldir = joinpath(pkgdir(Geom4hep), "examples");

world = read_gdml(joinpath(gdmldir, "polyhedra.gdml"))
[(d.volume.label, d.volume.shape) for d in world.daughters]

fig = Figure(size=(1000, 400))
draw!(LScene(fig[1, 1], show_axis=false), world)
draw!(LScene(fig[1, 2], show_axis=false), read_gdml(joinpath(gdmldir, "torus.gdml")))
draw!(LScene(fig[1, 3], show_axis=false), read_gdml(joinpath(gdmldir, "polycones.gdml")))
fig

cms = read_gdml(joinpath(gdmldir, "cms2018.gdml"), Float32)
[d.volume.label for d in cms[1].daughters][1:8]

barrel = cms[1, 7, 1]
world = world_volume(barrel)
fig = Figure(size=(800, 600))
draw!(LScene(fig[1, 1], show_axis=false), world; maxlevel=3)
fig

model = close_geometry(world)

fig = Figure(size=(1000, 500))
for axis in (3, 1)
    u, v, image, nsteps = xray(model, axis, 300)
    heatmap!(Axis(fig[1, axis == 3 ? 1 : 2], title="rays along $("xyz"[axis])", aspect=DataAspect()),
             u, v, image, colormap=:grayC)
end
fig
