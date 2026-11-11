#---Plane------------------------------------------------------------------------------------------
struct Plane{T<:AbstractFloat} <: AbstractShape{T}
    normal::Vector3{T}
    distance::T
end

#---Constructors-----------------------------------------------------------------------------------
function Plane{T}(normal::Vector3{T}, origin::Union{Point3{T},Vector3{T}}) where T<:AbstractFloat
    mag = norm(normal)
    Plane{T}(normal/mag, -dot(normal,origin)/mag)
end

function Plane{T}(θ, ϕ, origin::Point3{T}) where T<:AbstractFloat
    normal = Vector3{T}(sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    Plane{T}(normal, -dot(normal,origin))
end

function Plane{T}(p1::Point3{T}, p2::Point3{T}, p3::Point3{T}, p4::Point3{T}) where T<:AbstractFloat
    normal = (p4 - p1) × (p3 - p1)
    center = (p1 + p2 + p3 + p4)/4
    Plane{T}(normal, center)
end

#---Basic functionality----------------------------------------------------------------------------
function safety(plane::Plane{T}, point::Point3{T})::T where T<:AbstractFloat
    # Returns distance from point to plane
    # This is positive if the point is on the outside halfspace, negative otherwise.
    return dot(plane.normal, point) + plane.distance
end

function distance_to_in(plane::Plane{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    # The function returns a negative distance for points already inside or
    # direction going outwards (along the normal)
    d = T(-Inf)
    ndd = dot(normalize(dir), plane.normal)
    saf = safety(plane, point)
    ndd < 0 && saf > -tolerance(T) && (d = -saf/ndd)
    return d
end

function distance_to_out(plane::Plane{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    # The function returns infinity if the plane is not hit from inside, negative
    # if the point is outside
    d = T(Inf)
    ndd = dot(normalize(dir), plane.normal)
    saf = safety(plane, point)
    saf > tolerance(T) && (d = -T(Inf))
    ndd > 0 && saf < tolerance(T) && (d = -saf/ndd)
    return d
end

function z_limit(plane::Plane{T}, r::T, ϕ::T)::T where T<:AbstractFloat
    (; normal, distance) = plane
    return -(r * (cos(ϕ)*normal[1] + sin(ϕ)*normal[2]) + distance)/normal[3]
end

