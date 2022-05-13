using Geom4hep
using StaticArrays
using BenchmarkTools, Test, Printf
using GLMakie
using CUDA
using Profile

#-----------------------------------------------------------------------------------------------------------------------

function generateXRay(model::CuGeoModel, npoints::Number, view::Int=1; cuda::Bool=true) where T<:AbstractFloat
    world = model.volumes[1]
    lower, upper = extent(model.shapes[world.shapeIdx])
    dim = upper - lower
    if view == 1
        dim_a, dim_b = dim[2], dim[3]
    elseif view == 2
        dim_a, dim_b = dim[3], dim[1]
    elseif view == 3
        dim_a, dim_b = dim[1], dim[2]
    end
    pixel = round(sqrt(dim_a*dim_b/npoints), sigdigits=3)
    nx, ny = round(Int, dim_a/pixel), round(Int, dim_b/pixel)
    result = zeros(nx,ny)
    states = fill(CuNavigatorState{Float64}(1), nx, ny)
    if cuda && CUDA.functional()
        threads = (8,8)
        blocks = cld.((nx,ny),threads)
        cu_result = CuArray(result)
        cu_states = CuArray(states)
        cu_model = cu(model)
        CUDA.@sync @cuda threads=threads blocks=blocks k_generateXRay(cu_result, cu_states, cu_model, lower, pixel, view)
        return Array(cu_result)
    else
        generateXRay(result, states, model, lower, pixel, view)
        return result
    end
end

function  generateXRay(result::Matrix{T}, states::Matrix{CuNavigatorState{T}}, model::CuGeoModel, lower::Point3{T}, pixel::T, view::Int) where T<:AbstractFloat
    nx, ny = size(result)
    state = CuNavigatorState{T}(1)
    for i in 1:nx, j in 1:ny
        if view == 1
            point = Point3{T}(lower[1]+kTolerance(), lower[2]+(i-0.5)*pixel, lower[3]+(j-0.5)*pixel)
            dir = Vector3{T}(1,0,0)
        elseif view == 2
            point = Point3{T}(lower[1]+(j-0.5)*pixel, lower[2]+kTolerance(), lower[3]+(i-0.5)*pixel)
            dir = Vector3{T}(0,1,0)
        elseif view == 3
            point = Point3{T}(lower[1]+(i-0.5)*pixel, lower[2]+(j-0.5)*pixel, lower[3]+kTolerance())
            dir = Vector3{T}(0,0,1)
        end
        state = states[i,j]
        locateGlobalPoint!( model, state, point)
        mass::T =  0.0
        step::T = -1.
        while step != 0.0
            vol = model.volumes[state.currentVol]
            density = model.materials[vol.materialIdx].density
            step = computeStep!(model, state, point, dir, 1000.)
            point = point + dir * step
            mass += step * density
        end
        @inbounds result[i,j] = mass
    end
  
end
function k_generateXRay(result, states, model, lower::Point3{T}, pixel::T, view::Int) where T<:AbstractFloat
    nx, ny = size(result)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if i <= nx && j <= ny
        if view == 1
            point = Point3{T}(lower[1]+kTolerance(), lower[2]+(i-0.5)*pixel, lower[3]+(j-0.5)*pixel)
            dir = Vector3{T}(1,0,0)
        elseif view == 2
            point = Point3{T}(lower[1]+(j-0.5)*pixel, lower[2]+kTolerance(), lower[3]+(i-0.5)*pixel)
            dir = Vector3{T}(0,1,0)
        elseif view == 3
            point = Point3{T}(lower[1]+(i-0.5)*pixel, lower[2]+(j-0.5)*pixel, lower[3]+kTolerance())
            dir = Vector3{T}(0,0,1)
        end
        state = states[i,j]
        locateGlobalPoint!(model, state, point)
        mass::T =  0.0
        step::T = -1.0
        while step != 0.0
            vol = model.volumes[state.currentVol]
            density = model.materials[vol.materialIdx].density
            step = computeStep!(model, state, point, dir, 1000.)
            point = point + dir * step
            mass += step * density
        end
        @inbounds result[i,j] = mass
    end
    return
end

#-----build and generate image-----------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    world = processGDML("examples/boxes.gdml")
    model = fillCuGeometry(world)
    fig = Figure()
    @printf "generating x-projection\n"
    heatmap!(Axis(fig[1, 1], title = "X direction"), generateXRay(model, 1e4, 1), colormap=:grayC)
    @printf "generating y-projection\n"
    heatmap!(Axis(fig[2, 1], title = "Y direction"), generateXRay(model, 1e4, 2), colormap=:grayC)
    @printf "generating z-projection\n"
    heatmap!(Axis(fig[1, 2], title = "Z direction"), generateXRay(model, 1e4, 3), colormap=:grayC)
    draw(LScene(fig[2, 2]), world)
    save("xray-boxes.png", fig)
end
