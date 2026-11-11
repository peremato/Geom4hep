using Geom4hep
using StaticArrays

const idx = [2 3 1; 1 3 2; 1 2 3]
const rdx = [3 1 2; 1 3 2; 1 2 3]

# Generate a X-Ray of a geometry----------------------------------------------
function generate_image(vol::Volume{T}, npoints::Number, view::Int=1, d::T=T(0)) where T<:AbstractFloat
    world = world_volume(vol)
    lower, upper = extent(world.shape)
    ix, iy, iz = idx[view, :]
    xi, yi, zi = rdx[view, :]

    nx = ny = round(Int, sqrt(npoints))
    xrange = range(lower[ix], upper[ix], length = nx)
    yrange = range(lower[iy], upper[iy], length = ny)

    result = zeros(T, nx,ny)
    state = NavigatorState{T}(world)

    for (i,x) in enumerate(xrange), (j,y) in enumerate(yrange)
        _point = (x, y, d)
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
        locate_global_point!(state, point)
        result[i,j] = current_volume(state).material.density
    end
    return xrange, yrange, result
end

world = read_gdml("examples/trackML.gdml", Float64)
volume = world.daughters[2].volume.daughters[1].volume
#world = read_gdml("examples/polycones.gdml", Float64)
world = read_gdml("examples/cms2018.gdml", Float64)
volume = world

rx = generate_image(volume, 1e6, 1)
ry = generate_image(volume, 1e6, 2)
rz = generate_image(volume, 1e6, 3)
limits = (0, max(maximum(rx[3]), maximum(ry[3]), maximum(rz[3])))
using GLMakie
fig = Figure(size = (1000, 1000))
heatmap!(Axis(fig[1, 1], title = "X direction"), rx..., colormap=:grayC, colorrange=limits)
heatmap!(Axis(fig[1, 2], title = "Y direction"), ry..., colormap=:grayC, colorrange=limits)
heatmap!(Axis(fig[2, 1], title = "Z direction"), rz..., colormap=:grayC, colorrange=limits)
display(fig)