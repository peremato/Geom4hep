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




#=
<Float32>
==== Using shape: Box{Float32}(Float32[5.0, 10.0, 15.0]) with 10240 points and directions
CPU:   399.430 μs (0 allocations: 0 bytes)
GPU:   10.814 μs (3 allocations: 272 bytes)
==== Using shape: Tube{Float32}(5.0f0, 10.0f0, 10.0f0, 0.0f0, 6.2831855f0) with 10240 points and directions
CPU:   407.242 μs (0 allocations: 0 bytes)
GPU:   21.101 μs (3 allocations: 272 bytes)
==== Using shape: Cone{Float32}(5.0f0, 10.0f0, 7.0f0, 15.0f0, 10.0f0, 0.0f0, 6.2831855f0) with 10240 points and directions
CPU:   649.546 μs (0 allocations: 0 bytes)
GPU:   12.317 μs (3 allocations: 272 bytes)

<Float64>
==== Using shape: Box{Float64}([5.0, 10.0, 15.0]) with 10240 points and directions
CPU:   71.975 μs (0 allocations: 0 bytes)
GPU:   11.859 μs (3 allocations: 272 bytes)
==== Using shape: Tube{Float64}(5.0, 10.0, 10.0, 0.0, 6.283185307179586) with 10240 points and directions
CPU:   353.223 μs (0 allocations: 0 bytes)
GPU:   21.106 μs (3 allocations: 272 bytes)
==== Using shape: Cone{Float64}(5.0, 10.0, 7.0, 15.0, 10.0, 0.0, 6.283185307179586) with 10240 points and directions
CPU:   691.625 μs (0 allocations: 0 bytes)
GPU:   33.259 μs (3 allocations: 272 bytes)

=#


