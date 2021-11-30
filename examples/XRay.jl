using Geom4hep

# Geometry construction
function buildGeom(T::Type) 
    world = Volume{T}("World", Box{T}(100,100,100), Material("vacuum"; density=0.0))
    box1 = Volume{T}("box1", Box{T}(10,20,30), Material("iron"; density=7.0))
    box2 = Volume{T}("box2", Box{T}(4,4,4), Material("gold"; density=19.0))
    placeDaughter!(box1, Transformation3D{T}(0,0,0),box2)
    placeDaughter!(world, Transformation3D{T}( 50, 50, 50, RotXYZ{Float64}(π/4, 0, 0)), box1)
    placeDaughter!(world, Transformation3D{T}( 50, 50,-50, RotXYZ{Float64}(0, π/4, 0)), box1)
    placeDaughter!(world, Transformation3D{T}( 50,-50, 50, RotXYZ{Float64}(0, 0, π/4)), box1)
    placeDaughter!(world, Transformation3D{T}( 50,-50,-50, RotXYZ{Float64}(π/4, 0, 0)), box1)
    placeDaughter!(world, Transformation3D{T}(-50, 50, 50, RotXYZ{Float64}(0, π/4, 0)), box1)
    placeDaughter!(world, Transformation3D{T}(-50,-50, 50, RotXYZ{Float64}(0, 0, π/4)), box1)
    placeDaughter!(world, Transformation3D{T}(-50, 50,-50, RotXYZ{Float64}(π/4, 0, 0)), box1)
    placeDaughter!(world, Transformation3D{T}(-50,-50,-50, RotXYZ{Float64}(0, π/4, 0)), box1)
    trd1 = Volume{T}("trd1", Trd{T}(30,10,30,10,20), Material("water", density=10.0))
    placeDaughter!(world, Transformation3D{T}(0,0,0), trd1)
    return world
end

# Generate a X-Ray of a geometry
function generateXRay(world::Volume, npoints::Number, view::Symbol=:x)
    # Setup plane of points and results
    lower, upper = extent(world.shape)
    dim = upper - lower
    if view == :x
        shift = 0
    elseif view == :y
        shift = 1
    elseif view == :z
        shift = 2
    else
        error("Invalid projection: $view")
    end
    dummy, dim_a, dim_b = circshift(dim, shift)
    pixel = round(sqrt(dim_a*dim_b/npoints), sigdigits=3)
    nx, ny = round(Int, dim_a/pixel), round(Int, dim_b/pixel)
    result = zeros(nx,ny)

    state = NavigatorState{Float64}(world)
    dir = Vector3{Float64}(circshift([1,0,0], shift))
    point = Point3{Float64}(0,0,0)

    for i in 1:nx, j in 1:ny
        if shift == 0
            point = Point3{Float64}(lower[1]+kTolerance, lower[2]+(i-0.5)*pixel, lower[3]+(j-0.5)*pixel)
        elseif shift == 1
            point = Point3{Float64}(lower[1]+(j-0.5)*pixel, lower[2]+kTolerance, lower[3]+(i-0.5)*pixel)
        elseif shift == 2
            point = Point3{Float64}(lower[1]+(i-0.5)*pixel, lower[2]+(j-0.5)*pixel, lower[3]+kTolerance)
        end
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
            #move!(point, dir, step)
            point = point + dir * step
            mass += step * density
        end
        result[i,j] = mass   
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    #-----build and generate image
    world = buildGeom(Float64)
    image = generateXRay(world, 1e6)
    #-----Plot results
    using GLMakie
    heatmap(image, colormap=:grayC)
end

