#---Box-------------------------------------------------------------------------
struct Box{T<:AbstractFloat} <: AbstractShape{T}
    fDimensions::SVector{3,T} # the HALF lengths of the box
end

#---Contructors-----------------------------------------------------------------
Box(x::T,y::T,z::T) where T<:AbstractFloat = Box{T}(SVector{3,T}(x,y,z))
Box{T}(x::Number, y::Number, z::Number) where T<:AbstractFloat = Box{T}(SVector{3,T}(x,y,z))

#---Basic functions-------------------------------------------------------------
Base.getindex(b::Box, s::Symbol) = b.fDimensions[coordmap[s]]
Base.getindex(b::Box, i::Int64) = b.fDimensions[i]
#Base.getproperty(b::Box, s::Symbol) = (s in keys(coordmap) ? getfield(b, :fDimensions)[coordmap[s]] : getfield(b,s))
capacity(b::Box) = 8.0 * prod(b.fDimensions)
surface(b::Box) = 8.0(b.fDimensions[1] * b.fDimensions[2] + b.fDimensions[2] * b.fDimensions[3] + b.fDimensions[1] * b.fDimensions[3])
function extent(b::Box{T}) where T<:AbstractFloat
    (Point3{T}(-b.fDimensions), Point3{T}(b.fDimensions))
end

function normal(box::Box{T}, point::Point3{T}) where T
    safety = abs.(abs.(point) - box.fDimensions)
    safmin = minimum(safety)
    if safmin < tolerance(T)/2
        n = Vector3{T}(map( t -> t[1] - safmin < tolerance(T)/2 ? sign(t[2]) : 0.0, zip(safety, point)))
        normalize(n)
    else
        nothing    
    end
end

function inside(box::Box{T}, point::Point3{T})::Int64  where T<:AbstractFloat
    dist = -Inf
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > dist 
            dist = d 
        end
    end
    abs(dist) <= tolerance(T)/2 ? kSurface : dist < 0.0 ? kInside : kOutside
    #dist = maximum(abs.(point) - box.fDimensions)
    #isapprox(dist, 0.0, atol = tolerance(T)/2) ? kSurface : dist < 0.0 ? kInside : kOutside
end

function distance_to_out(box::Box{T}, point::Point3{T}, direction::Vector3{T})::T where T<:AbstractFloat
    safety = -Inf
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > safety 
            safety = d 
        end
    end
    if safety > tolerance(T)/2
        return -1.0
    end 
    dist = Inf
    for i in 1:3
        d = (copysign(box.fDimensions[i], direction[i]) - point[i])/direction[i]
        if d < dist
            dist = d
        end
    end
    return dist
    #safety = maximum(abs.(point) - box.fDimensions)
    #distance = minimum((copysign.(box.fDimensions, direction) - point) * (1.0 ./ direction))
    #safety > tolerance(T)/2 ? -1.0 : distance
end

function distance_to_in(box::Box{T}, point::Point3{T}, direction::Vector3{T},rcp_direction::Vector3{T}=inv.(direction))::T where T<:AbstractFloat
    
    max=Point3{T}(box.fDimensions)
    (distance,distout)= intersect_aabox_ray(-max,max,point,direction,rcp_direction)
    
    # The checks are as follows:
    # 1.) distance < -tolerance(T) would mean that we are already further in the volume than a micro step
    # 2.) distance >= distout Not sure if this is still needed
    # 3.) isnan(distance) NaN use to escape the intersect intersection routine and would cause problems
    # 4.) distout <= tolerance(T)/2) Checks to see if we are entering and immediatelly exiting
    (distance < -tolerance(T) || distance >= distout || isnan(distance) || distout <= tolerance(T)/2) ? Inf : distance
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
    v = 1/√3
    return Vector3{T}[(-v,-v,-v), ( v,-v,-v), (-v, v,-v), ( v, v,-v),
                      (-v,-v, v), ( v,-v, v), (-v, v, v), ( v, v, v)]
end
