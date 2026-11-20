using Geom4hep, Geom4hep.Units
using CairoMakie: CairoMakie, Figure, Axis, heatmap!, lines!, Colorbar

function build(T)
    world = Volume{T}("world", Box{T}(50cm, 50cm, 50cm), Material{T}("air"; density=0.0012))
    shell = Volume{T}("shell", Sphere{T}(36cm, 38cm), Material{T}("steel"; density=7.9))
    tank  = Volume{T}("tank", BooleanSubtraction(Box{T}(20cm, 20cm, 20cm), Tube{T}(0, 3cm, 25cm, 0, 2π)),
                      Material{T}("water"; density=1.0))
    ball  = Volume{T}("ball", Sphere{T}(0, 8cm), Material{T}("lead"; density=11.35))
    place_daughter!(tank, Transformation3D{T}(12cm, 0, 0), ball)
    place_daughter!(world, one(Transformation3D{T}), shell)
    place_daughter!(world, one(Transformation3D{T}), tank)
    return world
end
world = build(Float64);

model, index = close_geometry(world; with_index=true)
model

name(state) = in_world(state) ? index.volumes[state.volume].label : "outside"
[name(locate(model, Point3(x, 0.0, 0.0))) for x in (0.0, 12cm, 25cm, 37cm, 45cm, 60cm)]

function follow(model, point, dir)
    state = locate(model, point)
    crossings = []
    while in_world(state)
        step, next = compute_step(model, state, point, dir, Inf)
        push!(crossings, (volume=name(state), step=step))
        point += step * dir
        state = next
    end
    return crossings
end
follow(model, Point3(-49cm, 1cm, 0.0), Vector3(1.0, 0.0, 0.0))

follow(model, Point3(0.0, 0.0, -49cm), Vector3(0.0, 0.0, 1.0))

u, v, image, nsteps = xray(model, 3, 300)
fig = Figure(size=(900, 400))
ax = Axis(fig[1, 1], aspect=1, title="along z", xlabel="x (mm)", ylabel="y (mm)")
hm = heatmap!(ax, u, v, image, colormap=:grays)
u, v, image, nsteps = xray(model, 1, 300)
ax = Axis(fig[1, 2], aspect=1, title="along x", xlabel="y (mm)", ylabel="z (mm)")
heatmap!(ax, u, v, image, colormap=:grays)
Colorbar(fig[1, 3], hm, label="mass thickness (g/cm³ × mm)")
fig

model32 = close_geometry(build(Float32))
u, v, image32, _ = xray(model32, 1, 300)
maximum(abs, image32 - image) / maximum(image)
