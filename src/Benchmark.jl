#---Benchmarking functions-------------------------------------------------------------------------------------------------
using Printf
export benchmark_shape, fill_points!, fill_directions!

function fill_points!(points::Vector{Point3{T}}, shape::AbstractShape{T}, in_or_out::Integer=kInside) where T<:AbstractFloat
    low, hi = extent(shape)
    dim = hi - low
    if in_or_out == kOutside
        low -= dim/10  # increase the extent by -10% + 10%
        hi  += dim/10
        dim = hi - low
    end
    idx = 0
    while idx < length(points)
        point = low + rand(Vector3{T}) * dim
        if inside(shape, point) == in_or_out
            idx += 1
            points[idx] = point
        end
    end
end

function fill_directions!(dirs::Vector{Vector3{T}}, shape::AbstractShape{T}, points::Vector{Point3{T}}, bias::T=zero(T) ) where T<:AbstractFloat
    N = length(dirs)
    hits = 0
    for idx in 1:N
        dir = normalize(rand(Vector3{T}) + Vector3{T}(-.5,-.5,-.5))
        if distance_to_in(shape, points[idx], dir) < Inf
            hits += 1
        elseif hits/idx < bias
            while distance_to_in(shape, points[idx], dir) == Inf
                dir = normalize(rand(Vector3{T}) + Vector3{T}(-.5,-.5,-.5))
            end
            hits += 1
        end
        dirs[idx] = dir
    end
end


function time_function(func::Function, shape::AbstractShape{T}, points::Vector{Point3{T}}, dirs::Vector{Vector3{T}}) where T<:AbstractFloat
    N = length(points)
    # jit it eventually
    for i in 1:N
        func(shape, points[i], dirs[i])
    end
    start_t = time_ns()
    for i in 1:N
        func(shape, points[i], dirs[i])
    end
    end_t = time_ns()
    1e-9 * (end_t - start_t)
end
function time_function(func::Function, shape::AbstractShape{T}, points::Vector{Point3{T}}) where T<:AbstractFloat
    N = length(points)
    # jit it eventually
    for i in 1:N
        func(shape, points[i])
    end
    start_t = time_ns()
    for i in 1:N
        func(shape, points[i])
    end
    end_t = time_ns()
    1e-9 * (end_t - start_t)
end

function benchmark_shape(shape::AbstractShape{T}; npoints::Integer=1024, repetitions::Integer=1) where T<:AbstractFloat
    # Allocate Vectors of points and directions
    points = Vector{Point3{T}}(undef, npoints)
    dirs = Vector{Vector3{T}}(undef, npoints)
    bias = 0.8

    @printf "================================================================================\n"
    @printf "Running benchmark for %d points for %d repetitions and bias %f\n" npoints repetitions bias
    @printf "With shape: %s\n" "$shape"

    toin = 0.
    safin = 0.
    for r in 1:repetitions
        fill_points!(points, shape, kOutside)
        fill_directions!(dirs, shape, points, bias)
        toin += time_function(distance_to_in, shape, points, dirs)
        safin += time_function(safety_to_in, shape, points)
    end
    toout = 0.
    safout = 0.
    for r in 1:repetitions
        fill_points!(points, shape, kInside)
        fill_directions!(dirs, shape, points)
        toout += time_function(distance_to_out, shape, points, dirs)
        safout += time_function(safety_to_out, shape, points)
    end

    @printf " -- distance_to_in:  %9f s, safety_to_in:  %9f s, ratio: %6f\n" toin/repetitions safin/repetitions toin/safin
    @printf " -- distance_to_out: %9f s, safety_to_out: %9f s, ratio: %6f\n" toout/repetitions safout/repetitions toout/safout

end
