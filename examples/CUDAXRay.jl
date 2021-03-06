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

function generateXRay(model::CuGeoModel{Vector{Geom4hep.CuVolume{T}}}, npoints::Number, xview::Int=1; cuda::Bool=true) where T<:AbstractFloat
    world = model.volumes[1]
    lower, upper = extent(model.shapes[world.shapeIdx])
    idx = (idxs[xview,1], idxs[xview,2], idxs[xview,3])
    rdx = (rdxs[xview,1], rdxs[xview,2], rdxs[xview,3])
    ix, iy, iz = idx
    nx = ny = round(Int, sqrt(npoints))
    xrange = range(lower[ix], upper[ix], length = nx)
    yrange = range(lower[iy], upper[iy], length = ny)
    result = zeros(T,nx,ny)
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
    px = (upper[ix]-lower[ix])/nx
    py = (upper[iy]-lower[iy])/ny
    _dir = (0, 0, 1)
    dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi]) 
    state = CuNavigatorState{T}(1)
    for i in 1:nx, j in 1:ny
        _point = (lower[ix]+ (i - 0.5)*px, lower[iy]+ (j - 0.5)*py, lower[iz]+ kTolerance(T))
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
        locateGlobalPoint!( model, state, point)
        mass::T =  0.0
        step::T =  0.0
        while step >= 0.0
            vol = model.volumes[state.currentVol]
            density = model.materials[vol.materialIdx].density
            step = computeStep!(model, state, point, dir, T(1000))
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

    px = (upper[ix]-lower[ix])/(nx)
    py = (upper[iy]-lower[iy])/(ny)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if i <= nx && j <= ny
        _dir = (0, 0, 1)
        dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi])
        _point = (lower[ix]+ (i - T(0.5))*px, lower[iy]+ (j - T(0.5))*py, lower[iz]+ kTolerance(T))
        point = Point3{T}(_point[xi], _point[yi], _point[zi])

        state = CuNavigatorState{T}(1)
        locateGlobalPoint!(model, state, point)
        mass::T =  T(0)
        step::T =  T(0)
        while step >= 0
            vol = model.volumes[state.currentVol]
            density = model.materials[vol.materialIdx].density
            step = computeStep!(model, state, point, dir, T(1000))
            point = point + dir * step
            mass += step * density
        end
        @inbounds result[i,j] = mass
    end
    return
end

#-----build and generate image-----------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    #-----build detector---------------------------------------
    world = processGDML("examples/trackML.gdml", Float32)
    volume = world.daughters[2].volume.daughters[1].volume
    model = fillCuGeometry(getWorld(volume))
    #-----Generate maps----------------------------------------
    @printf "generating x-projection\n"
    rx = generateXRay(model, 1e6, 1)
    @printf "generating y-projection\n"
    ry = generateXRay(model, 1e6, 2)
    @printf "generating z-projection\n"
    rz = generateXRay(model, 1e6, 3)
    limits = (0, max(maximum(rx[3]), maximum(ry[3]), maximum(rz[3])))
    #----Plot the results--------------------------------------
    fig = Figure(resolution = (1000, 1000))
    @printf "ploting the results\n"
    heatmap!(Axis(fig[1, 1], title = "X direction"), rx..., colormap=:grayC, colorrange=limits)
    heatmap!(Axis(fig[2, 1], title = "Y direction"), ry..., colormap=:grayC, colorrange=limits)
    heatmap!(Axis(fig[1, 2], title = "Z direction"), rz..., colormap=:grayC, colorrange=limits)
    draw(LScene(fig[2, 2]), volume, 3)
    #display(fig)
    save("CuTrackML.png", fig)

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
