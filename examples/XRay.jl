using Geom4hep

# Generate a X-Ray of a geometry----------------------------------------------
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
            point = Point3{Float64}(lower[1]+kTolerance(), lower[2]+(i-0.5)*pixel, lower[3]+(j-0.5)*pixel)
        elseif shift == 1
            point = Point3{Float64}(lower[1]+(j-0.5)*pixel, lower[2]+kTolerance(), lower[3]+(i-0.5)*pixel)
        elseif shift == 2
            point = Point3{Float64}(lower[1]+(i-0.5)*pixel, lower[2]+(j-0.5)*pixel, lower[3]+kTolerance())
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
            point = point + dir * step
            mass += step * density
        end
        result[i,j] = mass   
    end
    return result
end

using GLMakie

if abspath(PROGRAM_FILE) == @__FILE__
    #-----build and generate image-----------------------------
    world = processGDML("examples/boxes.gdml")
    #-----Generate and Plot results----------------------------
    fig = Figure()
    heatmap!(LScene(fig[1, 1]), generateXRay(world, 1e6, :x), colormap=:grayC)
    heatmap!(LScene(fig[2, 1]), generateXRay(world, 1e6, :y), colormap=:grayC)
    heatmap!(LScene(fig[1, 2]), generateXRay(world, 1e6, :z), colormap=:grayC)
    draw(LScene(fig[2, 2]), world)
    display(fig)
end


