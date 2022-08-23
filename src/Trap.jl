#---Trap (General Trapezoid)------------------------------------------------------------------------
struct Trap{T<:AbstractFloat} <: AbstractShape{T}
    z::T    # Half-length along the z-axis
    θ::T    # Polar angle of the line joining the centres of the faces at -/+z
    ϕ::T    # Azimuthal angle of the line joing the centre of the face at -z to the centre of the face at +z
    y1::T   # Half-length along y of the face at -z
    x1::T   # Half-length along x of the side at y=-y1 of the face at -z
    x2::T   # Half-length along x of the side at y=+y1 of the face at -z
    α1::T   # Angle with respect to the y axis from the centre of the side at y=-y1 to the centre at y=+y1 of the face at -z
    y2::T   # Half-length along y of the face at +pDz
    x3::T   # Half-length along x of the side at y=-y2 of the face at +z
    x4::T   # Half-length along x of the side at y=+y2 of the face at +z
    α2::T   # Angle with respect to the y axis from the centre of the side at y=-y2 to the centre at y=+y2 of the face at +z

    coordinates::SVector{8,Point3{T}}       # Coordinates of 8 vertices 
    planes::SVector{4,Plane{T}}             # Planes of the 4 sides
end

#---Constructor------------------------------------------------------------------------------------
function Trap{T}(z, θ, ϕ, y1, x1, x2, α1, y2, x3, x4, α2) where T<:AbstractFloat
    ztanθcosϕ = z * tan(θ) * cos(ϕ)
    ztanθsinϕ = z * tan(θ) * sin(ϕ)
    y1tanα1 = y1 * tan(α1)
    y2tanα2 = y2 * tan(α2)
    pt = [Point3{T}(-ztanθcosϕ - y1tanα1 - x1, -ztanθsinϕ - y1, -z),
          Point3{T}(-ztanθcosϕ - y1tanα1 + x1, -ztanθsinϕ - y1, -z),
          Point3{T}(-ztanθcosϕ + y1tanα1 - x2, -ztanθsinϕ + y1, -z),
          Point3{T}(-ztanθcosϕ + y1tanα1 + x2, -ztanθsinϕ + y1, -z),
          Point3{T}( ztanθcosϕ - y2tanα2 - x3,  ztanθsinϕ - y2,  z),
          Point3{T}( ztanθcosϕ - y2tanα2 + x3,  ztanθsinϕ - y2,  z),
          Point3{T}( ztanθcosϕ + y2tanα2 - x4,  ztanθsinϕ + y2,  z),
          Point3{T}( ztanθcosϕ + y2tanα2 + x4,  ztanθsinϕ + y2,  z)]
    iface = [[1 5 6 2],[3 4 8 7],[1 3 7 5],[2 6 8 4]]
    planes = [Plane{T}(pt[f[1]], pt[f[2]], pt[f[3]], pt[f[4]]) for f in iface]
    Trap{T}(z, θ, ϕ, y1, x1, x2, α1, y2, x3, x4, α2, pt, planes)
end
function Trap{T}(x1, x2, y1, y2, z) where T<:AbstractFloat
    Trap{T}(z, 0, 0, y1, x1, x1, 0, y2, x2, x2, 0)
end

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, shape::Trap{T}) where T
    (; z, θ, ϕ, y1, x1, x2, α1, y2, x3, x4, α2) = shape
    print(io, "Trap{$T}",(z=z, θ=θ, ϕ=ϕ, y1=y1, x1=x1, x2=x2, α1=α1, y2=y2, x3=x3, x4=x4, α2=α2))
end

function GeometryBasics.coordinates(trap::Trap{T}, facets=6) where {T<:AbstractFloat}
    trap.coordinates
end

function GeometryBasics.faces(trap::Trap{T}, facets=6) where {T<:AbstractFloat}
    iface = ((1,5,6,2),(3,4,8,7),(1,3,7,5),(2,6,8,4),(1,2,4,3),(5,6,8,7))
    (QuadFace{Int64}(f...) for f in iface) 
end


#---Basic Functions---------------------------------------------------------------------------------
function capacity(trap::Trap{T}) where T<:AbstractFloat
    (; z, θ, ϕ, y1, x1, x2, α1, y2, x3, x4, α2) = trap
    z * ((x1 + x2 + x3 + x4) * (y1 + y2) + (x4 + x3 - x2 - x1) * (y2 - y1) / 3)
end

function surface(trap::Trap{T}) where T<:AbstractFloat
    p = trap.coordinates
    iface = ((1,5,6,2),(3,4,8,7),(1,3,7,5),(2,6,8,4),(1,2,4,3),(5,6,8,7))
    sum = T(0)
    for f in iface
        sum += norm((p[f[3]] - p[f[1]]) × (p[f[4]] - p[f[2]]))/2
    end
    return sum
end

function extent(trap::Trap{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    p = Point3{T}( max(trap.x1, trap.x2, trap.x3, trap.x4), max(trap.y1, trap.y2), trap.z)
    ( -p, p )
end

@override function inside(trap::Trap{T}, point::Point3{T})
    z = point[3]
    planes = trap.planes
    # inside z?
    outside = abs(z) > trap.z + kTolerance(T)/2
    cinside = abs(z) < trap.z - kTolerance(T)/2

    for i in 1:4
        s = safety(planes[i], point)  #  positive if the point is on the outside halfspace, negative otherwise.
        outside |= s > kTolerance(T)/2
        cinside &= s < -kTolerance(T)/2
    end
    cinside ? kInside : outside ? kOutside : kSurface
end

@override function safetyToOut(trap::Trap{T}, point::Point3{T}) where T<:AbstractFloat
    (; z, planes) = trap
    saf = z - abs(point[3])
    for i in 1:4
        s = -safety(planes[i], point)
        s < saf && (saf = s)
    end
    return saf
end

@override function safetyToIn(trap::Trap{T}, point::Point3{T}) where T<:AbstractFloat
    (; z, planes) = trap
    saf = abs(point[3]) - z
    for i in 1:4
        s = safety(planes[i], point)
        s > saf && (saf = s)
    end
    return saf
end

@override function distanceToOut(trap::Trap{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    (; z, planes) = trap
    dz = dir[3]

    # step 0: if point is outside any plane --> return -1, otherwise initialize at Infinity
    outside = abs(point[3]) > z + kTolerance(T)/2
    outside && return T(-1)
    distance = T(Inf)

    # Step 1: find range of distances along dir between Z-planes (smin, smax)
    dz != 0 && (distance = (copysign(z,dz) - point[3]) / dz)

    # Step 2: find distances for intersections with side planes.
    for i in 1:4
        dist = T(Inf)
        ndir = dot(dir, planes[i].normal)
        safe = safety(planes[i], point)

        safe > kTolerance(T) && return T(-1)
        ndir > 0 && safe < kTolerance(T) && (dist = -safe/ndir)
        dist < distance && (distance = dist)
    end
    return distance
end

@override function distanceToIn(trap::Trap{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    (; z, planes) = trap
    dz = dir[3]
    z = point[3]

    # step 1.a: input particle is moving away --> return infinity
    signZdir = copysign(1,dz)
    max = signZdir * trap.z - z
    signZdir * max < kTolerance(T)/2 && return Inf

    # step 1.b General case:
    smax = max/dz
    smin = -(signZdir * trap.z + z)/dz 

    # Step 2: find distances for intersections with side planes.
    done = false
    for i in 1:4
        ndir = dot(dir, planes[i].normal)
        safe = safety(planes[i], point)
        dist = -safe/ndir
        
        done |= (safe > kTolerance(T)/2 && ndir >= 0) || (safe > -kTolerance(T) && ndir > 0)

        posPoint = safe > -kTolerance(T)/2
        posDir   = ndir > 0
        !posPoint && posDir && dist < smax && (smax = dist)
        posPoint && !posDir && dist > smin && (smin = dist)
    end
    done && return T(Inf)
    distance = smin <= smax ? smin : Inf
    distance < -kTolerance(T)/2 && return T(-1)
    return distance
end 



