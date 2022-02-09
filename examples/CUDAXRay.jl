using Geom4hep
using CUDA
using BenchmarkTools, Test, Printf

function k_generateXRay(result, lower::Point3{T}, pixel::T, state::NavigatorState{T}, view::Symbol) where T<:AbstractFloat
    nx, ny = size(result)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if i <= nx && j <= ny
        if view == :x
            point = Point3{T}(lower[1]+kTolerance(), lower[2]+(i-0.5)*pixel, lower[3]+(j-0.5)*pixel)
            dir = Vector3{T}(1,0,0)
        elseif view == :y
            point = Point3{T}(lower[1]+(j-0.5)*pixel, lower[2]+kTolerance(), lower[3]+(i-0.5)*pixel)
            dir = Vector3{T}(0,1,0)
        elseif shift == :z
            point = Point3{T}(lower[1]+(i-0.5)*pixel, lower[2]+(j-0.5)*pixel, lower[3]+kTolerance())
            dir = Vector3{T}(0,0,1)
        end
        locateGlobalPoint!(state, point)
        mass::T =  0.0
        step::T = -1.0
        while step != 0.0
            density = currentVolume(state).material.density
            step = computeStep!(state, point, dir, 1000.)
            point = point + dir * step
            mass += step * density
        end
        @inbounds result[i,j] = mass
    end
    return
end

function generateXRay(world::Volume{T}, npoints::Number, view::Symbol=:x) where T<:AbstractFloat
    lower, upper = extent(world.shape)
    dim = upper - lower
    if view == :x
        dim_a, dim_b = dim[2], dim[3]
    elseif view == :y
        dim_a, dim_b = dim[3], dim[1]
    elseif view == :z
        dim_a, dim_b = dim[1], dim[2]
    end
    pixel = round(sqrt(dim_a*dim_b/npoints), sigdigits=3)
    nx, ny = round(Int, dim_a/pixel), round(Int, dim_b/pixel)
    state = NavigatorState{Float64}(world)
    result = CUDA.zeros(nx,ny)
    threads = (16,16)
    blocks = cld.((nx,ny),threads)
    CUDA.@sync @cuda threads=threads blocks=blocks k_generateXRay(result, lower, pixel, state, view)
    return Array(result)
end

#-----build and generate image-----------------------------
world = processGDML("examples/boxes.gdml")

