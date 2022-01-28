#---Benchmarking functions-------------------------------------------------------------------------------------------------
using Printf
export benchmarkShape

function fillPoints!(points::Vector{Point3{T}}, shape::AbstractShape{T}, in_or_out::Integer=kInside) where T<:AbstractFloat
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

function fillDirections!(dirs::Vector{Vector3{T}}, shape::AbstractShape{T}, points::Vector{Point3{T}}, bias::T=zero(T) ) where T<:AbstractFloat
    N = length(dirs)
    hits = 0
    for idx in 1:N
        dir = normalize(rand(Vector3{T}) + Vector3{T}(-.5,-.5,-.5))
        if distanceToIn(shape, points[idx], dir) < Inf
            hits += 1
        elseif hits/idx < bias
            while distanceToIn(shape, points[idx], dir) == Inf
                dir = normalize(rand(Vector3{T}) + Vector3{T}(-.5,-.5,-.5))
            end
            hits += 1
        end
        dirs[idx] = dir
    end
end


function timeFunction(func::Function, shape::AbstractShape{T}, points::Vector{Point3{T}}, dirs::Vector{Vector3{T}}) where T<:AbstractFloat
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
function timeFunction(func::Function, shape::AbstractShape{T}, points::Vector{Point3{T}}) where T<:AbstractFloat
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

function benchmarkShape(shape::AbstractShape{T}; Npoints::Integer=1024, Repetitions::Integer=1) where T<:AbstractFloat
    # Allocate Vectors of points and directions
    points = Vector{Point3{T}}(undef, Npoints)
    dirs = Vector{Vector3{T}}(undef, Npoints)
    bias = 0.8

    @printf "================================================================================\n"
    @printf "Running benchmark for %d points for %d repetitions and bias %f\n" Npoints Repetitions bias
    @printf "With shape: %s\n" "$shape"

    toin = 0.
    safin = 0.
    for r in 1:Repetitions
        fillPoints!(points, shape, kOutside)
        fillDirections!(dirs, shape, points, bias)
        toin += timeFunction(distanceToIn, shape, points, dirs)
        safin += timeFunction(safetyToIn, shape, points)
    end
    toout = 0.
    safout = 0.
    for r in 1:Repetitions
        fillPoints!(points, shape, kInside)
        fillDirections!(dirs, shape, points)
        toout += timeFunction(distanceToOut, shape, points, dirs)
        safout += timeFunction(safetyToOut, shape, points)
    end

    @printf " -- DistanceToIn:  %9f s, safetyToIn:  %9f s, DistanceToIn/SafetyToIn: %6f \n" toin/Repetitions safin/Repetitions toin/safin
    @printf " -- DistanceToOut: %9f s, safetyToOut: %9f s, DistanceToOut/SafetyToUot: %6f \n" toout/Repetitions safout/Repetitions toout/safout

end