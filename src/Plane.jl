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


function zLimit(plane::Plane{T}, r::T, ϕ::T)::T where T<:AbstractFloat
    (; normal, distance) = plane
    return -(r * (cos(ϕ)*normal[1] + sin(ϕ)*normal[2]) + distance)/normal[3]
end

