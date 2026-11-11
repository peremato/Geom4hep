#  Density maps of slices of a detector: each pixel shows the density of the material at a point of
#  a plane, found with `locate` on a flat GeometryModel.
#
#  Usage (from the repository root):
#     julia --project=examples -t auto examples/Render.jl

using Geom4hep, GLMakie
using Geom4hep: Point3        # Makie also exports a Point3

function density_slice(model::GeometryModel{T}, axis::Integer, position::Real, npixels::Integer; halfsize=Inf) where T
    lower, upper = shape_call(extent, model, model.volumes[1].shape)
    lower, upper = max.(lower, -T(halfsize)), min.(upper, T(halfsize))    # window of interest
    iu, iv = axis == 1 ? (2, 3) : axis == 2 ? (1, 3) : (1, 2)
    u = range(lower[iu], upper[iu], length = npixels)
    v = range(lower[iv], upper[iv], length = npixels)
    image = zeros(T, npixels, npixels)
    Threads.@threads for j in 1:npixels
        for i in 1:npixels
            c = zeros(T, 3); c[axis] = position; c[iu] = u[i]; c[iv] = v[j]
            state = locate(model, Point3{T}(c...))
            image[i, j] = in_world(state) ? model.materials[material_index(model, state)].density : 0
        end
    end
    return u, v, image
end

world = read_gdml(joinpath(@__DIR__, "cms2018.gdml"), Float64)
model = close_geometry(world)
println(model)

fig = Figure(size = (1500, 500))
for (axis, title) in ((1, "x = 0"), (2, "y = 0"), (3, "z = 0"))
    t = @elapsed u, v, img = density_slice(model, axis, 0, 1000; halfsize=10_000)   # ±10 m
    println("slice $title: $(round(t, digits=2)) s")
    heatmap!(Axis(fig[1, axis], title = title, aspect = DataAspect()), u, v, img, colormap = :grayC)
end
save("cms2018-slices.png", fig)
println("image saved in cms2018-slices.png")
