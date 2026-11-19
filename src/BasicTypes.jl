#--- basic stuff
"""
    AbstractShape{T<:AbstractFloat}

Supertype of all the shapes (solids). `T` is the floating point type of the shape parameters and of
the points and directions it works with (`Float64` or `Float32`). A shape is defined in its own
frame; it is positioned in a mother volume by the transformation of its placement.

Every shape implements [`inside`](@ref), [`distance_to_in`](@ref), [`distance_to_out`](@ref),
[`safety_to_in`](@ref), [`safety_to_out`](@ref), [`normal`](@ref), [`extent`](@ref),
[`capacity`](@ref) and [`surface`](@ref), and the mesh functions used for drawing
(`coordinates`, `faces`, [`mesh`](@ref)).
"""
abstract type AbstractShape{T<:AbstractFloat} end

"""
    AbstractMaterial{T<:AbstractFloat}

Supertype of [`Isotope`](@ref), [`Element`](@ref) and [`Material`](@ref).
"""
abstract type AbstractMaterial{T<:AbstractFloat} end

#---Shape functions (methods for each shape are in src/shapes)--------------------------------------
"""
    inside(shape, point) -> kInside | kSurface | kOutside

Classify `point` (in the frame of `shape`) as inside the shape, on its surface (within
[`tolerance`](@ref)) or outside it. The constants `kInside`, `kSurface` and `kOutside` are exported.
"""
function inside end

"""
    distance_to_in(shape, point, dir)

Distance along the unit direction `dir` from `point`, outside `shape`, to where the ray enters the
shape. `Inf` if the ray misses the shape. A point on the surface moving into the shape gets 0. The
result is not meaningful for points inside the shape (see [`inside`](@ref)).
"""
function distance_to_in end

"""
    distance_to_out(shape, point, dir)

Distance along the unit direction `dir` from `point`, inside `shape`, to where the ray leaves the
shape. A point on the surface moving out of the shape gets 0, and a point outside the shape gets a
negative value (-1).
"""
function distance_to_out end

"""
    safety_to_in(shape, point)

Lower bound of the distance from `point`, outside `shape`, to the shape (an isotropic safe step). It
is cheaper to compute than the exact distance; it is ≤ 0 for points inside the shape.
"""
function safety_to_in end

"""
    safety_to_out(shape, point)

Lower bound of the distance from `point`, inside `shape`, to its surface (an isotropic safe step). It
is ≤ 0 for points outside the shape.
"""
function safety_to_out end

"""
    extent(shape) -> (lower::Point3, upper::Point3)

Corners of an axis-aligned box, in the frame of the shape, that contains it.
"""
function extent end

"""
    capacity(shape)

Volume of the shape. For shapes without a closed formula it is estimated by sampling (see
[`estimate_capacity`](@ref)).
"""
function capacity end

"""
    surface(shape)

Area of the surface of the shape. For shapes without a closed formula it is estimated by sampling
(see [`estimate_surface`](@ref)).
"""
function surface end

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

"""
    @nounroll for/while ... end

Ask LLVM not to unroll the loop. On GPUs all calls are inlined, so unrolling a loop that calls large
functions multiplies the size of the kernel (see the CSG shapes).
"""
macro nounroll(ex)
    (ex isa Expr && ex.head in (:for, :while)) || error("@nounroll expects a for or while loop")
    push!(ex.args[2].args, Expr(:loopinfo, (Symbol("llvm.loop.unroll.disable"), 1)))
    return esc(ex)
end

#---constants
#const tolerance = 1e-9
const kInside  = 0
const kSurface = 1
const kOutside = 2
export kInside, kSurface, kOutside

"""
    tolerance(T=Float64)

Thickness of the surfaces of the shapes: a point closer than `tolerance(T)/2` to a surface is on
the surface ([`inside`](@ref) returns `kSurface`). 1 nm for `Float64` and 1 µm for `Float32`
(lengths in mm).
"""
tolerance(::Type{Float32}) = 1f-3
tolerance(::Type{Float64}) = 1e-9
tolerance() = tolerance(Float64)

"""
    push_tolerance(T=Float64)

Distance by which the navigator pushes a point beyond a boundary it has reached, to locate it in
the next volume.
"""
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
