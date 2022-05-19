using Revise
using Geom4hep
using Printf
using StaticArrays

const idx = [2 3 1; 1 3 2; 1 2 3]
const rdx = [3 1 2; 1 3 2; 1 2 3]

# Generate a X-Ray of a geometry----------------------------------------------
function generateXRay(world::Volume{T}, npoints::Number, view::Int=1) where T<:AbstractFloat
    # Setup plane of points and results
    lower, upper = extent(world.shape)
    ix, iy, iz = idx[view, :]
    xi, yi, zi = rdx[view, :]

    nx = ny = round(Int, sqrt(npoints))
    xrange = range(lower[ix], upper[ix], length = nx)
    yrange = range(lower[iy], upper[iy], length = ny)
    z = lower[iz]+kTolerance(T)
 
    result = zeros(nx,ny)

    state = NavigatorState{T}(world)
    _dir = (0, 0, 1)
    dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi]) 

    for (i,x) in enumerate(xrange), (j,y) in enumerate(yrange)
        _point = (x, y, z)
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
        locateGlobalPoint!(state, point)
        mass =  0.0
        step = -1.0
        while step != 0.0
            density = currentVolume(state).material.density
            step = computeStep!(state, point, dir, 1000.)
            if isnan(step)
                @show point, dir, state
                break
            end
            point = point + dir * step
            mass += step * density
        end
        result[i,j] = mass   
    end
    return xrange, yrange, result
end

using GLMakie

if abspath(PROGRAM_FILE) == @__FILE__
    #-----build detector---------------------------------------
    world = processGDML("examples/trackML.gdml")
    volume = world.daughters[2].volume.daughters[1].volume
    #-----Generate and Plot results----------------------------
    fig = Figure()
    @printf "generating x-projection\n"
    heatmap!(Axis(fig[1, 1], title = "X direction"), generateXRay(volume, 1e6, 1)..., colormap=:grayC)
    @printf "generating y-projection\n"
    heatmap!(Axis(fig[2, 1], title = "Y direction"), generateXRay(volume, 1e6, 2)..., colormap=:grayC)
    @printf "generating z-projection\n"
    heatmap!(Axis(fig[1, 2], title = "Z direction"), generateXRay(volume, 1e6, 3)..., colormap=:grayC)
    draw(LScene(fig[2, 2]), volume, 3)
    display(fig)
end


