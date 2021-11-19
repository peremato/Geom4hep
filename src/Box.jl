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
    if safmin < kTolerance/2
        n = Vec3{T}(map( t -> t[1] - safmin < kTolerance/2 ? sign(t[2]) : 0.0, zip(safety, point)))
        normalize(n)
    else
        nothing    
    end
end

function inside(box::Box{T}, point::Point3{T}) where T<:AbstractFloat
    dist = -Inf
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > dist 
            dist = d 
        end
    end
    abs(dist) <= kTolerance/2 ? kSurface : dist < 0.0 ? kInside : kOutside
    #dist = maximum(abs.(point) - box.fDimensions)
    #isapprox(dist, 0.0, atol = kTolerance/2) ? kSurface : dist < 0.0 ? kInside : kOutside
end

function distanceToOut(box::Box{T}, point::Point3{T}, direction::Vector3{T}) where T<:AbstractFloat
    safety = -Inf
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > safety 
            safety = d 
        end
    end
    if safety > kTolerance/2
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
    #safety > kTolerance/2 ? -1.0 : distance
end

function distanceToIn(box::Box{T}, point::Point3{T}, direction::Vector3{T}) where T<:AbstractFloat
    distsurf = Inf
    distance = -Inf
    distout = Inf
    for i in 1:3
        din  = (-copysign(box.fDimensions[i],direction[i]) - point[i])/direction[i]
        tout =   copysign(box.fDimensions[i],direction[i]) - point[i]
        dout = tout/direction[i]
        dsur = copysign(tout, direction[i])
        if din > distance 
            distance = din 
        end
        if dout < distout
            distout = dout
        end
        if dsur < distsurf
            distsurf = dsur
        end 
    end
    (distance >= distout || distout <= kTolerance/2 || abs(distsurf) <= kTolerance/2) ? Inf : distance
    #invdir = 1.0 ./ direction
    #tempIn  = -copysign.(box.fDimensions, direction) - point
    #tempOut =  copysign.(box.fDimensions, direction) - point
    #distance = maximum(tempIn * invdir)
    #distout  = minimum(tempOut * invdir)
    #distsurf = abs(minimum(copysign.(tempOut, direction)))
    #(distance >= distout || distout <= kTolerance/2 || distsurf <= kTolerance/2) ? Inf : distance
end

function safetyToOut(box::Box{T}, point::Point3{T}) where T<:AbstractFloat
    minimum(box.fDimensions - abs.(point))
end

function safetyToIn(box::Box{T}, point::Point3{T}) where T<:AbstractFloat
    maximum(abs.(point) - box.fDimensions)
end

function toMesh(box::Box{T}) where T<:AbstractFloat
    x, y, z = box.fDimensions
    positions = Point3{T}[(-x,-y,-z), ( x,-y,-z), (-x, y,-z), ( x, y,-z),
                          (-x,-y, z), ( x,-y, z), (-x, y, z), ( x, y, z)]
    faces = TriangleFace{Int64}[(1,2,4), (4,3,1), (1,3,7), (7,5,1), (1,5,6), (6,2,1), 
                                (8,6,5), (5,7,8), (8,7,3), (3,4,8), (8,4,2), (2,6,8)] 
    return GeometryBasics.Mesh(positions, faces)
end

#-------------------------------------------------------------------------------
#---TBox------------------------------------------------------------------------
const box_faces = [(1,2,4), (4,3,1), (1,3,7), (7,5,1), (1,5,6), (6,2,1), 
                   (8,6,5), (5,7,8), (8,7,3), (3,4,8), (8,4,2), (2,6,8)]

struct TBox{T<:AbstractFloat} <: AbstractShape{T}
    dx::T # the HALF lengths of the box
    dy::T
    dz::T
    triangles::Vector{Triangle{T}}
    function TBox{T}(dx::Number, dy::Number, dz::Number) where T<:AbstractFloat
        vertices = Point3{T}[(-dx,-dy,-dz), ( dx,-dy,-dz), (-dx, dy,-dz), ( dx, dy,-dz),
                             (-dx,-dy, dz), ( dx,-dy, dz), (-dx, dy, dz), ( dx, dy, dz)]
        triangles = [Triangle{T}(vertices[i],vertices[j], vertices[k]) for (i,j,k) in box_faces]
        new(dx, dy, dz, triangles)
    end
end

function Base.show(io::IO, box::TBox{T}) where T<:AbstractFloat
    print(io, "TBox{$T}",(box.dx, box.dy, box.dz))
end

#=
function distanceToOut(box::TBox{T}, point::Point3{T}, direction::Vector3{T}) where T<:AbstractFloat
    for triangle in box.triangles
        dist, ok = intersect(point, direction, triangle)
        ok && return dist
    end
    -1.0
end
=#

function extent(b::TBox{T}) where T<:AbstractFloat
    (Point3{T}(-b.dx, -b.dy, -b.dz), Point3{T}(b.x, b.y, b.z))
end

