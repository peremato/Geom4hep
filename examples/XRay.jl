using Revise
using Geom4hep
using Printf
using StaticArrays
using TimerOutputs

const idx = [2 3 1; 1 3 2; 1 2 3]
const rdx = [3 1 2; 1 3 2; 1 2 3]

# Generate a X-Ray of a geometry----------------------------------------------
function generateXRay(nav::AbstractNavigator, vol::Volume{T}, npoints::Number, view::Int=1) where T<:AbstractFloat
    world = getWorld(vol)
    # Setup plane of points and results
    lower, upper = extent(world.shape)
    ix, iy, iz = idx[view, :]
    xi, yi, zi = rdx[view, :]

    #nx = ny = round(Int, sqrt(npoints))
    nx = npoints
    ny = trunc(Int, nx * (upper[iy] - lower[iy])/(upper[ix] - lower[ix]))
    xrange = range(lower[ix], upper[ix], length = nx+1)[1:nx]
    yrange = range(lower[iy], upper[iy], length = ny+1)[1:ny]
    z = lower[iz] + kPushTolerance(T)
 
    result = zeros(T, nx,ny)
    state = NavigatorState{T}(world, nav)
    _dir = (0, 0, 1)
    dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi]) 

    for (i,x) in enumerate(xrange), (j,y) in enumerate(yrange)
        _point = (x, y, z)
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
        #@timeit Geom4hep.to "point" locateGlobalPoint!(state, point)
        locateGlobalPoint!(state, point)
        mass = T(0)
        step = T(0)
        nsteps = 0
        maxsteps = 1000
        while isInVolume(state) && nsteps < maxsteps
            density = currentVolume(state).material.density
            #@timeit Geom4hep.to "computeStep" step = computeStep!(state, point, dir, T(10000))
            step = computeStep!(state, point, dir, T(10000))
            if step == -1
                error("step == -1 may indicate a navigation error. \nstate = $state \npoint = $point \ndir =   $dir" )
            end
            if step > 0
                point = point + dir * step
                mass += step * density
            end
            nsteps += 1
            if nsteps == maxsteps
                error("Reached max number of steps. \nstep = $step \nstate = $state \npoint = $point \ndir =   $dir" )
            end
        end
        result[i,j] = -log(nsteps+1)
    end
    return xrange, yrange, result
end

using GLMakie

#if abspath(PROGRAM_FILE) == @__FILE__
    #-----build detector---------------------------------------
    #full = processGDML("examples/trackML.gdml", Float64)
    #volume = full[2,1]
 
    full = processGDML("examples/cms2018.gdml", Float64)
    volume = full[1,7,1]
 
    world = getWorld(volume)
    nav = BVHNavigator(world)
    #nav = TrivialNavigator(world)

    #-----Generate and Plot results----------------------------
    @printf "generating x-projection\n"
    rx = generateXRay(nav, world, 1000, 1)
    @printf "generating y-projection\n"
    ry = generateXRay(nav, world, 1000, 2)
    @printf "generating z-projection\n"
    rz = generateXRay(nav, world, 1000, 3)
    limits = (min(minimum(rx[3]), minimum(ry[3]), minimum(rz[3])), max(maximum(rx[3]), maximum(ry[3]), maximum(rz[3])))
 
    #----Plot the results--------------------------------------
    fig = Figure(resolution = (1000, 1000))
    @printf "ploting the results\n"
    heatmap!(Axis(fig[1, 1], title = "X direction"), rx..., colormap=:grayC, colorrange=limits)
    heatmap!(Axis(fig[2, 1], title = "Y direction"), ry..., colormap=:grayC, colorrange=limits)
    heatmap!(Axis(fig[1, 2], title = "Z direction"), rz..., colormap=:grayC, colorrange=limits)
    #draw!(LScene(fig[2, 2]), volume; wireframe=true, maxlevel=2)
    #display(fig)
    save("cms2018.png", fig)
#end
