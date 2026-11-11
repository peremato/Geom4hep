#---Points and vectors------------------------------------------------------------------------------
#   Geom4hep owns its point/vector types, so that the affine arithmetic below does not change the
#   behaviour of other packages (e.g. GeometryBasics or Makie). They are isbits FieldVectors, so all
#   StaticArrays functionality (indexing, broadcasting, dot, cross, norm, ...) is available and they
#   can be stored in GPU arrays.
#---------------------------------------------------------------------------------------------------
"""
    Point3{T}(x, y, z)

A position in 3D space. `Point3 - Point3` gives a [`Vector3`](@ref), `Point3 ± Vector3` a `Point3`.
"""
struct Point3{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

"""
    Vector3{T}(x, y, z)

A displacement or direction in 3D space. Transformations rotate but do not translate vectors.
"""
struct Vector3{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

"""
    Point2{T}(x, y)

A position in 2D space (used for computations in the transverse plane).
"""
struct Point2{T} <: FieldVector{2,T}
    x::T
    y::T
end

const Vector2 = SVector{2}

#---Results of generic StaticArrays operations keep the Geom4hep type--------------------------------
StaticArrays.similar_type(::Type{<:Point3}, ::Type{T}, ::Size{(3,)}) where T = Point3{T}
StaticArrays.similar_type(::Type{<:Vector3}, ::Type{T}, ::Size{(3,)}) where T = Vector3{T}
StaticArrays.similar_type(::Type{<:Point2}, ::Type{T}, ::Size{(2,)}) where T = Point2{T}

#---Affine arithmetic--------------------------------------------------------------------------------
@inline Base.:-(a::Point3, b::Point3) = Vector3(a.x - b.x, a.y - b.y, a.z - b.z)
@inline Base.:+(a::Point3, b::Vector3) = Point3(a.x + b.x, a.y + b.y, a.z + b.z)
@inline Base.:+(a::Vector3, b::Point3) = b + a
@inline Base.:-(a::Point3, b::Vector3) = Point3(a.x - b.x, a.y - b.y, a.z - b.z)
# Sum of two points (only meaningful when followed by a scaling, e.g. a mid point (a + b)/2)
@inline Base.:+(a::Point3, b::Point3) = Vector3(a.x + b.x, a.y + b.y, a.z + b.z)

#---Element-wise products (scaling of extents, random points in a box, ...)--------------------------
@inline Base.:*(a::Vector3, b::Vector3) = Vector3(a.x * b.x, a.y * b.y, a.z * b.z)
@inline Base.:*(a::Point3, b::Vector3) = Point3(a.x * b.x, a.y * b.y, a.z * b.z)
@inline Base.:*(a::Vector3, b::Point3) = b * a

#---dot and cross with explicit evaluation order (bit-for-bit identical results to the original
#   implementation; some shape algorithms are sensitive to exact cancellations)
@inline LinearAlgebra.dot(a::Vector3, b::Vector3) = a.x * b.x + a.y * b.y + a.z * b.z
@inline LinearAlgebra.cross(a::Vector3, b::Vector3) = Vector3(a.y * b.z - a.z * b.y, -a.x * b.z + a.z * b.x, a.x * b.y - a.y * b.x)

unitize(v::Vector3) = v/sqrt(v⋅v)

#---Conversions to/from GeometryBasics points (for meshing and drawing)------------------------------
_gbpoint(p::StaticVector{3,T}) where T = GeometryBasics.Point3{T}(p[1], p[2], p[3])
