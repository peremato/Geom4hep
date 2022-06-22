#--- basic stuff
abstract type AbstractShape{T<:AbstractFloat} end
abstract type AbstractMaterial{T<:AbstractFloat} end

function GeometryBasics.Tesselation(shape::AbstractShape{T}, nvertices::NTuple{N,<:Integer}) where {T,N}
    return Tesselation{3,T,typeof(shape),N}(shape, Int.(nvertices))
end
  
const Vector3 = SVector{3}
const Vector2 = SVector{2}
Base.:*(x::Vector3, y::Vector3) = Vector3(x[1]*y[1], x[2]*y[2], x[3]*y[3])
Base.:*(x::Point3, y::Vector3) = Point3(x[1]*y[1], x[2]*y[2], x[3]*y[3])
Base.:*(x::Vector3, n::Number) = Vector3(x[1]*n, x[2]*n, x[3]*n)
Base.:*(n::Number, x::Vector3) = x * n
Base.:/(p::Point3, n::Number) = Point3(p[1]/n, p[2]/n, p[3]/n)
LinearAlgebra.:⋅(x::Vector3, y::Vector3) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
Base.:-(p1::Point3, p2::Point3) = Vector3(p1[1]-p2[1], p1[2]-p2[2], p1[3]-p2[3])
Base.:+(p1::Point3, p2::Point3) = Vector3(p1[1]+p2[1], p1[2]+p2[2], p1[3]+p2[3])
LinearAlgebra.:×(x::Vector3, y::Vector3) = Vector3(x[2]*y[3]-x[3]*y[2], -x[1]*y[3]+x[3]*y[1], x[1]*y[2]-x[2]*y[1])
unitize(v::Vector3) = v/(v⋅v)

#---constants
#const kTolerance = 1e-9
const kInside  = 0
const kSurface = 1
const kOutside = 2
export kInside, kSurface, kOutside

kTolerance(::Type{Float32}) = 1f-3
kTolerance(::Type{Float64}) = 1e-9
kTolerance() = kTolerance(Float64)
nonzero(x::T) where T<:AbstractFloat = x == zero(T) ? eps(T) : x
