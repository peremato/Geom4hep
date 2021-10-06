using Geom4hep

# Geometry construction
function buildGeom(T::Type)
    world = Volume("world", Box{T}(100,100,100), Material("vacuum", 0.0))
    box1 = Volume("box1", Box{T}(10,20,30), Material("iron", 7.0))
    box2 = Volume("box2", Box{T}(4,4,4), Material("gold", 19.0))
    placeDaughter!(box1, Transformation3D{T}(0,0,0),box2)
    placeDaughter!(world, Transformation3D{T}(50,50,50, RotXYZ{Float64}(0, π/4, 0)), box1)
    placeDaughter!(world, Transformation3D{T}(-50,-50,-50, RotXYZ{Float64}(π/4, 0, 0)), box1)
    return world
end

# Generate a X-Ray of a geometry
function generateXRay(world::Volume, npoints)
    # Setup plane of points and results
    lower, upper = extent(world.shape)
    dim = upper - lower
    pixel = round(sqrt(dim[2]*dim[3]/npoints), sigdigits=3)
    nx, ny = round(Int, dim[2]/pixel), round(Int, dim[3]/pixel)
    result = zeros(nx,ny)

    state = NavigatorState(world)
    dir = Vector3{Float64}(1,0,0)

    for i in 1:nx, j in 1:ny
        point = Point3{Float64}(lower[1]+kTolerance, lower[2]+(i-0.5)*pixel, lower[3]+(j-0.5)*pixel)
        locateGlobalPoint!(state, world, point)
        mass =  0.0
        step = -1.0
        while step != 0.0
            density = currentVolume(state).material.density
            step = computeStep!(state, point, dir, 1000.)
            if isnan(step)
                @show point, dir, state
                break
            end
            point += step * dir
            mass += step * density
        end
        result[i,j] = mass   
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    #-----build and generate image
    world = buildGeom(Float64)
    @time image = generateXRay(world, 1e6)
    #-----Plot results
    using GLMakie
    heatmap(image, colormap=:grayC)
end

