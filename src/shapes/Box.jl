#---Box-------------------------------------------------------------------------
struct Box{T<:AbstractFloat} <: AbstractShape{T}
    fDimensions::SVector{3,T} # the HALF lengths of the box
end

#---Contructors-----------------------------------------------------------------
Box(x::T,y::T,z::T) where T<:AbstractFloat = Box{T}(SVector{3,T}(x,y,z))
Box{T}(x::Number, y::Number, z::Number) where T<:AbstractFloat = Box{T}(SVector{3,T}(x,y,z))

#---Basic functions-------------------------------------------------------------
Base.getindex(b::Box, i::Integer) = b.fDimensions[i]
capacity(b::Box) = 8 * prod(b.fDimensions)
surface(b::Box) = 8 * (b.fDimensions[1] * b.fDimensions[2] + b.fDimensions[2] * b.fDimensions[3] + b.fDimensions[1] * b.fDimensions[3])
function extent(b::Box{T}) where T<:AbstractFloat
    (Point3{T}(-b.fDimensions), Point3{T}(b.fDimensions))
end

function normal(box::Box{T}, point::Point3{T}) where T<:AbstractFloat
    safety = abs.(abs.(point) - box.fDimensions)
    safmin = minimum(safety)
    safmin < tolerance(T)/2 || return zero(Vector3{T})
    n = Vector3{T}(ntuple(i -> safety[i] - safmin < tolerance(T)/2 ? sign(point[i]) : zero(T), 3))
    return normalize(n)
end

function inside(box::Box{T}, point::Point3{T})::Int64  where T<:AbstractFloat
    dist = T(-Inf)
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > dist 
            dist = d 
        end
    end
    abs(dist) <= tolerance(T)/2 ? kSurface : dist < 0 ? kInside : kOutside
end

function distance_to_out(box::Box{T}, point::Point3{T}, direction::Vector3{T})::T where T<:AbstractFloat
    safety = T(-Inf)
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > safety 
            safety = d 
        end
    end
    if safety > tolerance(T)/2
        return T(-1)
    end 
    dist = T(Inf)
    for i in 1:3
        d = (copysign(box.fDimensions[i], direction[i]) - point[i])/direction[i]
        if d < dist
            dist = d
        end
    end
    return dist
end

function distance_to_in(box::Box{T}, point::Point3{T}, direction::Vector3{T},rcp_direction::Vector3{T}=inv.(direction))::T where T<:AbstractFloat
    
    # Checks by face, on the distance to the face (not along the ray, which is large at grazing angles):
    # outside or on a face and moving away from it: no entry; further inside than the surface: Inf
    deepest = typemin(T)
    for i in 1:3
        safe = abs(point[i]) - box.fDimensions[i]
        pd = point[i] * direction[i]
        ((safe > tolerance(T)/2 && pd >= 0) || (safe > -tolerance(T)/2 && pd > 0)) && return T(Inf)
        deepest = max(deepest, safe)
    end
    deepest < -tolerance(T) && return T(Inf)

    bmax = Point3{T}(box.fDimensions)
    (distance,distout)= intersect_aabox_ray(-bmax, bmax, point, direction, rcp_direction)
    # A point on the surface moving in gets distance <= 0. The other checks:
    # - distance >= distout or NaN: the ray misses the box
    # - distout <= tolerance(T)/2: entering and immediately exiting
    (distance >= distout || isnan(distance) || distout <= tolerance(T)/2) ? T(Inf) : distance
end

function safety_to_out(box::Box{T}, point::Point3{T}) where T<:AbstractFloat
    minimum(box.fDimensions - abs.(point))
end

function safety_to_in(box::Box{T}, point::Point3{T}) where T<:AbstractFloat
    maximum(abs.(point) - box.fDimensions)
end

function GeometryBasics.coordinates(box::Box{T}, facets=6) where {T<:AbstractFloat}
    x, y, z = box.fDimensions
    return Point3{T}[(-x,-y,-z), ( x,-y,-z), (-x, y,-z), ( x, y,-z),
                     (-x,-y, z), ( x,-y, z), (-x, y, z), ( x, y, z)]
end

function GeometryBasics.faces(box::Box{T}, facets=6) where {T<:AbstractFloat}
    iface = ((1,5,6,2),(3,4,8,7),(1,3,7,5),(2,6,8,4),(1,2,4,3),(5,6,8,7))
    (QuadFace{Int64}(f...) for f in iface)
end

function GeometryBasics.normals(box::Box{T}, facets=24) where {T<:AbstractFloat}
    v = 1/√T(3)
    return Vector3{T}[(-v,-v,-v), ( v,-v,-v), (-v, v,-v), ( v, v,-v),
                      (-v,-v, v), ( v,-v, v), (-v, v, v), ( v, v, v)]
end
