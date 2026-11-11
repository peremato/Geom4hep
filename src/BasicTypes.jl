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

"""
    normal(shape, point)

Unit outward normal of `shape` at `point`. On edges and corners the normals of the surfaces meeting
there are averaged. Returns the zero vector if `point` is not on the surface (within tolerance).
"""
function normal end

#---constants
#const tolerance = 1e-9
const kInside  = 0
const kSurface = 1
const kOutside = 2
export kInside, kSurface, kOutside

tolerance(::Type{Float32}) = 1f-3
tolerance(::Type{Float64}) = 1e-9
tolerance() = tolerance(Float64)

push_tolerance(::Type{Float32}) = 1f-3
push_tolerance(::Type{Float64}) = 1e-6
push_tolerance() = push_tolerance(Float64)

nonzero(x::T) where T<:AbstractFloat = x == zero(T) ? eps(T) : x

#---Sampling estimates of capacity and surface (used for shapes without a closed formula)------------
#   Points come from the R_d low-discrepancy sequence (M. Roberts, 2018): deterministic results and a
#   faster convergence than pseudo-random sampling.
@inline function _rd_point(n::Integer, ::Val{D}) where D
    φ = D == 3 ? 1.2207440845519 : D == 5 ? 1.1347241384015194 : error("dimension $D not supported")
    ntuple(i -> mod(0.5 + n / φ^i, 1.0), D)
end

"""
    estimate_capacity(shape; npoints=1_000_000)

Volume of `shape` estimated from the fraction of `npoints` sample points of its bounding box that are
inside it. Typical accuracy with the default number of points is better than 0.1%.
"""
function estimate_capacity(shape::AbstractShape{T}; npoints::Integer=1_000_000) where T
    lo, hi = extent(shape)
    d = hi - lo
    nin = 0
    for n in 1:npoints
        u = _rd_point(n, Val(3))
        p = Point3{T}(lo[1] + u[1] * d[1], lo[2] + u[2] * d[2], lo[3] + u[3] * d[3])
        inside(shape, p) != kOutside && (nin += 1)
    end
    return T(prod(Float64.(d)) * nin / npoints)
end

"""
    estimate_surface(shape; npoints=1_000_000)

Surface area of `shape` estimated by sampling segments of length δ with isotropic directions in its
(enlarged) bounding box of volume V: the probability that a segment crosses the surface is A δ / 2V.
Typical accuracy with the default number of points is about 1% (slightly underestimated).
"""
function estimate_surface(shape::AbstractShape{T}; npoints::Integer=1_000_000) where T
    lo, hi = extent(shape)
    δ = minimum(hi - lo) / 100
    lo, hi = lo - Vector3{T}(δ, δ, δ), hi + Vector3{T}(δ, δ, δ)
    d = hi - lo
    crossings = 0
    for n in 1:npoints
        u = _rd_point(n, Val(5))
        p = Point3{T}(lo[1] + u[1] * d[1], lo[2] + u[2] * d[2], lo[3] + u[3] * d[3])
        cosθ = 2u[4] - 1
        sinθ = sqrt(1 - cosθ^2)
        ϕ = 2π * u[5]
        q = p + δ * Vector3{T}(sinθ * cos(ϕ), sinθ * sin(ϕ), cosθ)
        (inside(shape, p) == kOutside) != (inside(shape, q) == kOutside) && (crossings += 1)
    end
    return T(2 * prod(Float64.(d)) * crossings / (npoints * δ))
end
