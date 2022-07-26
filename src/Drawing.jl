using Makie
using Colors
using Printf

colors = colormap("Grays", 8)

#---Draw a Volume---------------------------------------------------------------
function draw(s::LScene, vol::Volume{T}, t::Transformation3D{T}, level::Int64, maxlevel::Int64) where T
    m = GeometryBasics.mesh(Tesselation(vol.shape, 64), facetype=typeof(first(GeometryBasics.faces(vol.shape))))
    if isone(t)
        mesh!(s, m, color=colors[level], transparency=true, ambient=0.7, visible = level == 1 ? false : true)
    else
        points = GeometryBasics.coordinates(m)
        faces  = GeometryBasics.faces(m)
        map!(c -> c * t, points, points)
        mesh!(s, points, faces, color=colors[level], transparency=true, ambient=0.7)
    end
    if level < maxlevel
        for daughter in vol.daughters
            draw(s, daughter.volume, t * daughter.transformation, level+1, maxlevel)
        end
    end
end

function draw(s::LScene, vol::Volume{T}, maxlevel::Int64=999) where T
    draw(s, vol, one(Transformation3D{T}), 1, maxlevel)
    display(s)
end

function draw(vol::Volume, maxlevel::Int64=999)
    fig = Figure()
    s = LScene(fig[1,1])
    draw(s, vol, one(Transformation3D{Float64}), 1, maxlevel)
    display(fig)
end

#---Draw a Shape---------------------------------------------------------------
function draw!(s::LScene, shape::AbstractShape; wireframe::Bool=false)
    m = GeometryBasics.mesh(Tesselation(shape, 32), facetype=typeof(first(faces(shape))))
    if wireframe
        wireframe!(s, m)
    else
        mesh!(s, m, color=:gray, transparency=true, ambient=0.7)
    end
    return s
end

function draw(shape::AbstractShape; wireframe::Bool=false)
    fig = Figure()
    s = LScene(fig[1,1])
    draw!(s, shape; wireframe)
    display(fig)
end

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

