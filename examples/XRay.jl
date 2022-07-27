using Revise
using Geom4hep
using Printf
using StaticArrays

const idx = [2 3 1; 1 3 2; 1 2 3]
const rdx = [3 1 2; 1 3 2; 1 2 3]

# Generate a X-Ray of a geometry----------------------------------------------
function generateXRay(vol::Volume{T}, npoints::Number, view::Int=1) where T<:AbstractFloat
    world = getWorld(vol)
    # Setup plane of points and results
    lower, upper = extent(world.shape)
    ix, iy, iz = idx[view, :]
    xi, yi, zi = rdx[view, :]

    nx = ny = round(Int, sqrt(npoints))
    xrange = range(lower[ix], upper[ix], length = nx)
    yrange = range(lower[iy], upper[iy], length = ny)
    z = lower[iz] + kTolerance(T)
 
    result = zeros(T, nx,ny)

    state = NavigatorState{T}(world)
    _dir = (0, 0, 1)
    dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi]) 

    for (i,x) in enumerate(xrange), (j,y) in enumerate(yrange)
        _point = (x, y, z)
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
        locateGlobalPoint!(state, point)
        mass = T(0)
        step = T(0)
        nsteps = 0
        while isInVolume(state) && nsteps < 100
            density = currentVolume(state).material.density
            step = computeStep!(state, point, dir, T(1000))
            if step > 0
                point = point + dir * step
                mass += step * density
            end
            nsteps += 1
        end
        result[i,j] = mass
    end
    return xrange, yrange, result
end

using GLMakie

if abspath(PROGRAM_FILE) == @__FILE__
    #-----build detector---------------------------------------
    world = processGDML("examples/trackML.gdml", Float64)
    volume = world.daughters[2].volume.daughters[1].volume
    #-----Generate and Plot results----------------------------
    @printf "generating x-projection\n"
    rx = generateXRay(volume, 1e6, 1)
    @printf "generating y-projection\n"
    ry = generateXRay(volume, 1e6, 2)
    @printf "generating z-projection\n"
    rz = generateXRay(volume, 1e6, 3)
    #limits = (0, max(maximum(rx[3]), maximum(ry[3]), maximum(rz[3])))
    #----Plot the results--------------------------------------
    fig = Figure(resolution = (1000, 1000))
    @printf "ploting the results\n"
    heatmap!(Axis(fig[1, 1], title = "X direction"), rx..., colormap=:grayC, colorrange=(0,2.))
    heatmap!(Axis(fig[2, 1], title = "Y direction"), ry..., colormap=:grayC, colorrange=(0,2.))
    heatmap!(Axis(fig[1, 2], title = "Z direction"), rz..., colormap=:grayC, colorrange=(0,maximum(rz[3])))
    draw(LScene(fig[2, 2]), volume, 3)
    #display(fig)
    save("trackML.png", fig)
end

