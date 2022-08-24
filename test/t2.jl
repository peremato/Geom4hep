using Revise
using GeometryBasics
using Geom4hep
using LinearAlgebra
using BenchmarkTools
using GLMakie

tube = Tube{Float64}(1,2,2,0,π)
#tube = Tube{Float64}(0,50,50,0,2π)
ctest10 = Cone{Float64}(20, 60, 80, 140, 100, 10deg, 300deg)

function drawDistanceToOut(shape::AbstractShape{T}, N::Integer) where T<:AbstractFloat
    low, hi = extent(shape)
    dim = hi - low
    points = (low + rp * dim for rp in rand(Vector3{Float64}, N))
    result = Vector{Point3{Float64}}()
    for point in points
        if inside(shape, point) == kInside
            dir = rand(Vector3) + Vector3(-.5,-.5,-.5)
            push!(result, point + dir * distanceToOut(shape, point, dir))
        end
    end
    fig = Figure()
    s = LScene(fig[1, 1])
    scatter!(s, result, color=:black, makersize=20)
    scatter!(s, [low, hi], color=:blue, markersize=50)
    display(fig)
    return s
end

function drawDistanceToIn(shape::AbstractShape{T}, N::Integer) where T<:AbstractFloat
    low, hi = extent(shape)
    dim = hi - low
    low -= dim/10.
    hi  += dim/10.
    dim = hi - low
    points = (low + rp * dim for rp in rand(Vector3{Float64}, N))
    dirs = normalize.(rand(Vector3{Float64}, N) .+ Ref(Vector3(-.5,-.5,-.5)))
    result = Vector{Point3{Float64}}()
    for (point, dir) in zip(points, dirs)
        if inside(shape, point) == kOutside
            dist = distanceToIn(shape, point, dir)
            if dist != Inf
                push!(result, point + dir * dist)
            end
        end
    end
    fig = Figure()
    s = LScene(fig[1, 1])
    scatter!(s, result, color=:black, markerspace=SceneSpace, markersize=dim[1]/500)
    scatter!(s, [low, hi], color=:blue, markerspace=SceneSpace, markersize=dim[1]/100)
    display(fig)
    return s
end






