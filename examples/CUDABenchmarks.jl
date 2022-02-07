using Geom4hep
using CUDA
using BenchmarkTools, Test, Printf


const N = 1024*10

#---CPU implementation----------------------------------------------
function cpu_d2out!(dist, shape, points, dirs)
    N = length(points)
    for i in 1:N
        @inbounds dist[i] = distanceToOut(shape, points[i], dirs[i])
    end
    return
end

#---GPU implementation----------------------------------------------
function gpu_d2out!(dist, shape, points, dirs)
    N = length(points)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:N
        @inbounds dist[i] = distanceToOut(shape, points[i], dirs[i])
    end
    return
end
function bench_gpu!(dist, shape, points, dirs)
    numblocks = ceil(Int, length(points)/256)
    CUDA.@sync begin
        @cuda threads=256 blocks=numblocks gpu_d2out!(dist, shape, points, dirs)
    end
end

for T in (Float32, Float64)

    cone = Cone{T}(5, 10, 7, 15, 10, 0, 2π)
    tube = Tube{T}(5,10, 10, 0, 2π)
    box  = Box{T}(5, 10, 15)

    points = Vector{Point3{T}}(undef, N)
    dirs = Vector{Vector3{T}}(undef, N)

    for shape in (box, tube, cone)
        fillPoints!(points, shape, kInside)
        fillDirections!(dirs, shape, points)
        dist = fill(zero(T), N)
        cpu_d2out!(dist, shape, points, dirs)

        dist_d = CUDA.fill(zero(T), N)
        points_d = CuArray(points)
        dirs_d = CuArray(dirs)

        @printf "==== Using shape: %s with %d points and directions\n" "$shape" N

        numblocks = ceil(Int, N/256)
        @cuda threads=256 blocks=numblocks gpu_d2out!(dist_d, shape, points_d, dirs_d)
        @test dist ≈ Array(dist_d)
    
        @printf "CPU: "
        @btime cpu_d2out!($dist, $shape, $points, $dirs)
        @printf "GPU: "
        @btime bench_gpu!($dist_d, $shape, $points_d, $dirs_d)
    end

end
