using Makie
using Colors
using Printf

import GeometryBasics: coordinates, faces, mesh, Mesh

colors = colormap("Grays", 16)

global drawn = 0

GeometryBasics.mesh(shape::AbstractShape{T}) where T = GeometryBasics.mesh(Tesselation(shape, 64), facetype=typeof(first(GeometryBasics.faces(shape))))

#---Draw a Volume---------------------------------------------------------------
function draw!(s::LScene, vol::Volume{T}, t::Transformation3D{T}, level::Int64, wireframe::Bool, bvh::Bool, maxlevel::Int64) where T
    m = GeometryBasics.mesh(vol.shape)
    if ! isone(t)
        points = GeometryBasics.coordinates(m)
        faces  = GeometryBasics.faces(m)
        map!(c -> c * t, points, points)
        m = GeometryBasics.Mesh(points, faces)
    end
    if wireframe
        wireframe!(s, m, color=colors[16-level], visible = level == 1 ? false : true)
    else
        mesh!(s, m, color=colors[level], transparency=true, ambient=0.7, visible = level == 1 ? false : true)
    end
    global drawn
    drawn += 1
    if drawn % 100 == 0; @show drawn; end
    if level < maxlevel
        for daughter in vol.daughters
            if level == 1 && bvh
                draw!(s, AABB(daughter))
            end
            draw!(s, daughter.volume, daughter.transformation * t, level+1, wireframe, bvh, maxlevel) 
        end
        if level == 1 && bvh
            draw!(s, buildBVH(vol.daughters))
        end
    end
end

function draw!(s::LScene, vol::Volume{T}; wireframe::Bool=false,  maxlevel::Int64=999) where T
    draw!(s, vol, one(Transformation3D{T}), 1, wireframe, false, maxlevel)
    display(s)
end

function draw(vol::Volume; wireframe::Bool=false, bvh::Bool=false, maxlevel::Int64=999)
    fig = Figure(resolution=(1280, 720))
    s = LScene(fig[1,1])
    draw!(s, vol, one(Transformation3D{Float64}), 1, wireframe, bvh, maxlevel)
    display(fig)
end

#---Draw a Shape---------------------------------------------------------------
function draw!(s::LScene, shape::AbstractShape; wireframe::Bool=false)
    m = GeometryBasics.mesh(shape)
    if wireframe
        wireframe!(s, m)
    else
        mesh!(s, m, color=colors[2], transparency=true, ambient=0.7)
    end
    return s
end

function draw(shape::AbstractShape; wireframe::Bool=false)
    fig = Figure()
    s = LScene(fig[1,1])
    draw!(s, shape; wireframe)
    display(fig)
end

#---Draw a AABB------------------------------------------------------------------
function draw!(s::LScene, aabb::AABB; color::Symbol=:blue)
    wireframe!(s, Mesh(coordinates(aabb), faces(aabb)), color=color)
end

function draw!(s::LScene, bvh::BVH; color::Symbol=:green)
    if isempty(bvh.children)
        draw!(s, bvh.aabb, color=color)
    else
        for c in bvh.children
            draw!(s, c, color=color)
        end
    end
end

function draw(bvh::BVH; color::Symbol=:green)
    fig = Figure()
    s = LScene(fig[1,1])
    draw!(s, bvh; color=color)
    display(fig)
end

#---Testing functions------------------------------------------------------------

function drawDistanceToOut(shape::AbstractShape{T}, N::Integer) where T<:AbstractFloat
    low, hi = extent(shape)
    dim = hi - low
    points = (low + rp * dim for rp in rand(Vector3{Float64}, N))
    result = Vector{Point3{Float64}}()
    for point in points
        if inside(shape, point) == kInside
            dir = normalize(rand(Vector3) + Vector3(-.5,-.5,-.5))
            push!(result, point + dir * distanceToOut(shape, point, dir))
        end
    end
    fig = Figure()
    s = LScene(fig[1, 1])
    scatter!(s, result, color=:black, markersize=1)
    scatter!(s, [low, hi], color=:blue, markersize=10)
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
    scatter!(s, result, color=:black, markersize=1)
    scatter!(s, [low, hi], color=:blue, markersize=10)
    display(fig)
    return s
end

