#---MultiUnion--------------------------------------------------------------------------------------
#   Union of any number of shapes, each one with its own placement (equivalent to G4MultiUnion).
#   Long chains of binary unions are converted to a MultiUnion when a flat GeometryModel is built,
#   which keeps the nesting of boolean shapes shallow. The storage of the components is a type
#   parameter: Vectors for shapes built by the user, views into a flat GeometryModel otherwise.
struct MultiUnion{T<:AbstractFloat, VS<:AbstractVector, VT<:AbstractVector{Transformation3D{T}}, VA<:AbstractVector{AABB{T}}} <: AbstractShape{T}
    shapes::VS              # component shapes
    transformations::VT     # placement of each component (`transformation * point` gives the component frame)
    aabbs::VA               # bounding box of each component in the frame of the MultiUnion (with tolerance)
end

#---Constructors-------------------------------------------------------------------------------------
function MultiUnion{T}(shapes::AbstractVector, transformations::AbstractVector{Transformation3D{T}}) where T<:AbstractFloat
    @assert length(shapes) == length(transformations) "Different number of shapes and transformations"
    tol = Vector3{T}(tolerance(T), tolerance(T), tolerance(T))
    aabbs = [(e = transform_extent(extent(s), t); AABB{T}(e[1] - tol, e[2] + tol)) for (s, t) in zip(shapes, transformations)]
    MultiUnion{T,Vector{TreeOperand{T}},Vector{Transformation3D{T}},Vector{AABB{T}}}(collect(TreeOperand{T}, shapes), collect(transformations), aabbs)
end
MultiUnion{T}() where T<:AbstractFloat = MultiUnion{T}(TreeOperand{T}[], Transformation3D{T}[])
MultiUnion(shapes::AbstractVector{<:AbstractShape{T}}, transformations::AbstractVector{Transformation3D{T}}) where T = MultiUnion{T}(shapes, transformations)

"""
    add_component!(mu::MultiUnion, shape, transformation=one(Transformation3D{T}))

Add a component to a MultiUnion built by the user.
"""
function add_component!(mu::MultiUnion{T}, shape::AbstractShape{T}, transformation::Transformation3D{T}=one(Transformation3D{T})) where T
    tol = Vector3{T}(tolerance(T), tolerance(T), tolerance(T))
    e = transform_extent(extent(shape), transformation)
    push!(mu.shapes, shape)
    push!(mu.transformations, transformation)
    push!(mu.aabbs, AABB{T}(e[1] - tol, e[2] + tol))
    return mu
end

ncomponents(mu::MultiUnion) = length(mu.shapes)

#---Printing and Plotting----------------------------------------------------------------------------
function Base.show(io::IO, mu::MultiUnion{T}) where T
    print(io, "MultiUnion{$T} with ", ncomponents(mu), " components")
end

function GeometryBasics.mesh(mu::MultiUnion{T}; facets::Integer=64) where T
    meshes = [(m = mesh(mu.shapes[i]; facets); Mesh(map(c -> _gbpoint(Point3{T}(c) * mu.transformations[i]), coordinates(m)), faces(m)))
              for i in 1:ncomponents(mu)]
    merge(meshes)
end

#---Basic functions----------------------------------------------------------------------------------
function extent(mu::MultiUnion{T})::Tuple{Point3{T},Point3{T}} where T
    tol = Vector3{T}(tolerance(T), tolerance(T), tolerance(T))
    lo = Point3{T}(typemax(T), typemax(T), typemax(T))
    hi = Point3{T}(typemin(T), typemin(T), typemin(T))
    for bb in mu.aabbs
        lo = min.(lo, bb.min + tol)
        hi = max.(hi, bb.max - tol)
    end
    return (lo, hi)
end

capacity(mu::MultiUnion; npoints::Integer=1_000_000) = estimate_capacity(mu; npoints)
surface(mu::MultiUnion; npoints::Integer=1_000_000) = estimate_surface(mu; npoints)

function inside(mu::MultiUnion{T}, point::Point3{T})::Int64 where T
    first_surface = 0
    for i in 1:ncomponents(mu)
        inside(mu.aabbs[i], point) || continue
        t = mu.transformations[i]
        pos = inside(mu.shapes[i], t * point)
        pos == kInside && return kInside
        if pos == kSurface
            if first_surface == 0
                first_surface = i
            else
                # two components touching with opposite normals: the point is inside the union  -)(-
                t1 = mu.transformations[first_surface]
                n1 = normal(mu.shapes[first_surface], t1 * point) * t1
                n2 = normal(mu.shapes[i], t * point) * t
                n1 ⋅ n2 < 0 && return kInside
            end
        end
    end
    return first_surface == 0 ? kOutside : kSurface
end

function normal(mu::MultiUnion{T}, point::Point3{T})::Vector3{T} where T
    # normal of a component on whose surface the point is, unless the point is inside another one
    n = zero(Vector3{T})
    for i in 1:ncomponents(mu)
        inside(mu.aabbs[i], point) || continue
        t = mu.transformations[i]
        pos = inside(mu.shapes[i], t * point)
        pos == kInside && return zero(Vector3{T})
        pos == kSurface && n ⋅ n == 0 && (n = normal(mu.shapes[i], t * point) * t)
    end
    return n
end

function distance_to_in(mu::MultiUnion{T}, point::Point3{T}, dir::Vector3{T})::T where T
    invdir = inv.(dir)
    distance = T(Inf)
    for i in 1:ncomponents(mu)
        intersect(mu.aabbs[i], point, dir, invdir) || continue
        t = mu.transformations[i]
        d = distance_to_in(mu.shapes[i], t * point, t * dir)
        d < distance && (distance = d)
    end
    return distance
end

function distance_to_out(mu::MultiUnion{T}, point::Point3{T}, dir::Vector3{T})::T where T
    # Move to the furthest exit point of the components containing the current point, until the
    # point is outside all of them (overlapping or touching components are crossed). After the
    # first step the probe point is pushed beyond the exit point; with a grazing crossing it can
    # still be within the tolerance of a surface, and then the direction decides.
    distance = zero(T)        # distance to the last exit point found
    offset = zero(T)          # distance to the probe point (distances are measured from there)
    npoint = point
    for iteration in 1:ncomponents(mu) + 1
        dmax = T(-1)
        for i in 1:ncomponents(mu)
            inside(mu.aabbs[i], npoint) || continue
            t = mu.transformations[i]
            lpoint = t * npoint
            pos = inside(mu.shapes[i], lpoint)
            pos == kOutside && continue
            ldir = t * dir
            # a probe point on the surface of a component is only followed if the ray enters it
            iteration > 1 && pos == kSurface && normal(mu.shapes[i], lpoint) ⋅ ldir >= 0 && continue
            d = distance_to_out(mu.shapes[i], lpoint, ldir)
            d > dmax && (dmax = d)
        end
        if dmax < 0                        # outside all components
            return iteration == 1 ? T(-1) : distance
        end
        distance = offset + dmax
        dmax == T(Inf) && return distance
        offset = distance + push_tolerance(T)
        npoint = point + offset * dir
    end
    return distance
end

function safety_to_in(mu::MultiUnion{T}, point::Point3{T})::T where T
    # distance to the union is the smallest distance to the components (all are lower bounds)
    safety = T(Inf)
    for i in 1:ncomponents(mu)
        s = safety_to_in(mu.shapes[i], mu.transformations[i] * point)
        s < safety && (safety = s)
    end
    return safety
end

function safety_to_out(mu::MultiUnion{T}, point::Point3{T})::T where T
    # a ball inside one component is inside the union
    safety = zero(T)
    for i in 1:ncomponents(mu)
        s = safety_to_out(mu.shapes[i], mu.transformations[i] * point)
        s > safety && (safety = s)
    end
    return safety
end

#---Operand types for shapes built by the user (union-split dispatch)--------------------------------
const TreeOperand{T} = Union{BaseShape{T}, BooleanShape{T}, MultiUnion{T}}
