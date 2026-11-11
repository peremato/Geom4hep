#--- basic stuff
abstract type AbstractShape{T<:AbstractFloat} end
abstract type AbstractMaterial{T<:AbstractFloat} end

"""
    mesh(shape::AbstractShape; facets=64)

Surface mesh of a shape (`GeometryBasics.Mesh`), e.g. for visualization. `facets` controls the
number of segments used to approximate curved surfaces.
"""
function GeometryBasics.mesh(shape::AbstractShape{T}; facets::Integer=64) where T
    points = [_gbpoint(p) for p in GeometryBasics.coordinates(shape, facets)]
    GeometryBasics.Mesh(points, collect(GeometryBasics.faces(shape, facets)))
end

#---constants
#const kTolerance = 1e-9
const kInside  = 0
const kSurface = 1
const kOutside = 2
export kInside, kSurface, kOutside

kTolerance(::Type{Float32}) = 1f-3
kTolerance(::Type{Float64}) = 1e-9
kTolerance() = kTolerance(Float64)

kPushTolerance(::Type{Float32}) = 1f-3
kPushTolerance(::Type{Float64}) = 1e-6
kPushTolerance() = kPushTolerance(Float64)

nonzero(x::T) where T<:AbstractFloat = x == zero(T) ? eps(T) : x
