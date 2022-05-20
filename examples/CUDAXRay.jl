using Geom4hep
using StaticArrays
using BenchmarkTools, Test, Printf
if haskey(ENV, "DISPLAY")
    using GLMakie
else
    using CairoMakie
end
using CUDA
using Profile

#-----------------------------------------------------------------------------------------------------------------------

const idxs =  [2 3 1; 1 3 2; 1 2 3]
const rdxs =  [3 1 2; 1 3 2; 1 2 3]

function generateXRay(model::CuGeoModel, npoints::Number, view::Int=1; cuda::Bool=true) where T<:AbstractFloat
    world = model.volumes[1]
    lower, upper = extent(model.shapes[world.shapeIdx])
    idx = (idxs[view,1], idxs[view,2], idxs[view,3])
    rdx = (rdxs[view,1], rdxs[view,2], rdxs[view,3])
    ix, iy, iz = idx
    nx = ny = round(Int, sqrt(npoints))
    xrange = range(lower[ix], upper[ix], length = nx)
    yrange = range(lower[iy], upper[iy], length = ny)
    result = zeros(nx,ny)
    if cuda && CUDA.functional()
        threads = (8,8)
        blocks = cld.((nx,ny),threads)
        cu_result = CuArray(result)
        cu_model = cu(model)
        CUDA.@sync @cuda threads=threads blocks=blocks k_generateXRay(cu_result, cu_model, lower, upper, idx, rdx)
        return xrange, yrange, Array(cu_result)
    else
        generateXRay(result, model, lower, upper, idx, rdx)
        return xrange, yrange, result
    end
end

function  generateXRay(result::Matrix{T}, model::CuGeoModel, lower::Point3{T}, upper::Point3{T}, idx::Tuple, rdx::Tuple) where T<:AbstractFloat
    nx, ny = size(result)
    ix, iy, iz = idx
    xi, yi, zi = rdx
    px = (upper[ix]-lower[ix])/(nx-1)
    py = (upper[iy]-lower[iy])/(ny-1)
    _dir = (0, 0, 1)
    dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi]) 
    state = CuNavigatorState{T}(1)
    for i in 1:nx, j in 1:ny
        _point = (lower[ix]+ i*px, lower[iy]+ j*py, lower[iz]+ kTolerance(T))
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
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
function k_generateXRay(result, model, lower::Point3{T}, upper::Point3{T}, idx::Tuple, rdx::Tuple) where T<:AbstractFloat
    nx, ny = size(result)

    ix, iy, iz = idx
    xi, yi, zi = rdx

    px = (upper[ix]-lower[ix])/(nx-1)
    py = (upper[iy]-lower[iy])/(ny-1)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if i <= nx && j <= ny
        _dir = (0, 0, 1)
        dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi])

        _point = (lower[ix]+ i*px, lower[iy]+ j*py, lower[iz]+ kTolerance(T))
        point = Point3{T}(_point[xi], _point[yi], _point[zi])

        state = CuNavigatorState{T}(1)
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
    world = processGDML("examples/trackML.gdml")
    topvol = world.daughters[2].volume.daughters[1].volume
    model = fillCuGeometry(topvol)
    fig = Figure()
    @printf "generating x-projection\n"
    heatmap!(Axis(fig[1, 1], title = "X direction"), generateXRay(model, 1e6, 1)..., colormap=:grayC)
    @printf "generating y-projection\n"
    heatmap!(Axis(fig[2, 1], title = "Y direction"), generateXRay(model, 1e6, 2)..., colormap=:grayC)
    @printf "generating z-projection\n"
    heatmap!(Axis(fig[1, 2], title = "Z direction"), generateXRay(model, 1e6, 3)..., colormap=:grayC)
    draw(LScene(fig[2, 2]), topvol, 3)
    save("xray-boxes.png", fig)
    if CUDA.functional()
        t1 = @benchmark generateXRay(model, 1e6, 3, cuda=true) seconds=60
        @show t1
        t2 = @benchmark generateXRay(model, 1e6, 3, cuda=false) seconds=60
        @show t2
        @printf "speedup [CPU/GPU]= %f\n" mean(t2).time/mean(t1).time
    else
        t2 = @benchmark generateXRay(model, 1e6, 3, cuda=false) seconds=60
        @show t2
    end
end
