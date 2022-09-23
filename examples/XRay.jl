using Revise
using Geom4hep
using Printf
using StaticArrays

const idx = [2 3 1; 1 3 2; 1 2 3]
const rdx = [3 1 2; 1 3 2; 1 2 3]

# Generate a X-Ray of a geometry----------------------------------------------
function generateXRay(nav::AbstractNavigator, vol::Volume{T}, npoints::Number, view::Int=1) where T<:AbstractFloat
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
    state = NavigatorState{T}(world, nav)
    _dir = (0, 0, 1)
    dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi]) 

    for (i,x) in enumerate(xrange), (j,y) in enumerate(yrange)
        _point = (x, y, z)
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
        locateGlobalPoint!(state, point)
        mass = T(0)
        step = T(0)
        nsteps = 0
        maxsteps = 1000
        while isInVolume(state) && nsteps < maxsteps
            density = currentVolume(state).material.density
            step = computeStep!(state, point, dir, T(1000))
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
        result[i,j] = mass
    end
    return xrange, yrange, result
end

const full = processGDML("examples/trackML.gdml", Float64)
const volume = full[2,1]
const world = getWorld(volume)
const nav = BVHNavigator(world)

const full2 = processGDML("examples/cms2018.gdml", Float64)
const volume2 = full2[1,7,1]
const world2 = getWorld(volume2)
const nav2 = BVHNavigator(world2)

function mainbench()
    @time generateXRay(nav, world, 1e4, 1);
    @time generateXRay(nav2, world2, 1e4, 1);
end
