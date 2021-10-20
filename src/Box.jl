#---Box-------------------------------------------------------------------------
struct Box{T<:AbstractFloat} <: AbstractShape{T}
    fDimensions::SVector{3,T} # the HALF lengths of the box
end
Box(x::T,y::T,z::T) where T<:AbstractFloat = Box{T}(SVector{3,T}(x,y,z))
Box{T}(x::Number, y::Number, z::Number) where T<:AbstractFloat = Box{T}(SVector{3,T}(x,y,z))
Base.getindex(b::Box, s::Symbol) = b.fDimensions[coordmap[s]]
Base.getindex(b::Box, i::Int64) = b.fDimensions[i]
#Base.getproperty(b::Box, s::Symbol) = (s in keys(coordmap) ? getfield(b, :fDimensions)[coordmap[s]] : getfield(b,s))
capacity(b::Box) = 8.0 * prod(b.fDimensions)
surface(b::Box) = 8.0(b.fDimensions[1] * b.fDimensions[2] + b.fDimensions[2] * b.fDimensions[3] + b.fDimensions[1] * b.fDimensions[3])
extent(b::Box{T}) where T<:AbstractFloat = (Point3{T}(-b.fDimensions), Point3{T}(b.fDimensions))

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
    orig, dest = extent(box)
    GeometryBasics.mesh(Rect3{T}(orig, dest-orig))
end