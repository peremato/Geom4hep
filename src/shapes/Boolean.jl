#---Basic Shape--------------------------------------------------------------------------------------
const BaseShape{T<:AbstractFloat} = Union{NoShape{T},
                                          Box{T},
                                          Trd{T},
                                          Trap{T},
                                          Tube{T},
                                          Cone{T},
                                          Polycone{T},
                                          Polyhedra{T},
                                          Torus{T},
                                          CutTube{T},
                                          Sphere{T}}

#---Boolean Shapes----------------------------------------------------------------------------------
#   The types of the operands are parameters (L, R). Shapes built by the user use a small Union of
#   shape types (see TreeOperand), while a flat GeometryModel uses references into its shape arrays.
struct BooleanShape{T<:AbstractFloat, OP, L, R} <: AbstractShape{T}
    left::L                                      # the mother (or left) volume A in unplaced form
    right::R                                     # (or right) volume B in placed form, acting on A with a boolean operation
    transformation::Transformation3D{T}          # placement of "right" with respect of "left" 
    left_aabb::AABB{T}
    right_aabb::AABB{T}
end

"""
    BooleanUnion(left, right, placement=one(Transformation3D{T}))
    BooleanSubtraction(left, right, placement=one(Transformation3D{T}))
    BooleanIntersection(left, right, placement=one(Transformation3D{T}))

Boolean combination of two shapes (as G4UnionSolid, G4SubtractionSolid and G4IntersectionSolid):
the points in `left` or `right`, in `left` and not in `right`, or in both. `right` is placed in the
frame of `left` with `placement` (as in [`place_daughter!`](@ref)). The operands can themselves be
boolean shapes.
"""
const BooleanUnion{T} = BooleanShape{T, :union}
"See [`BooleanUnion`](@ref)."
const BooleanSubtraction{T} =  BooleanShape{T, :subtraction}
"See [`BooleanUnion`](@ref)."
const BooleanIntersection{T} = BooleanShape{T, :intersection}

# Since all the other constructors only use 3 arguments we can sneak the AABB calc in here.
# The boxes are enlarged by the tolerance, so that the points on the surface of an operand that touches
# its bounding box are not discarded by the quick checks.
function BooleanShape{T,OP}(left, right, place) where {T, OP}
    tol = Vector3{T}(tolerance(T), tolerance(T), tolerance(T))
    lo, hi = extent(left)
    left_aabb = AABB{T}(lo - tol, hi + tol)
    lo, hi = transform_extent(extent(right), place)
    right_aabb = AABB{T}(lo - tol, hi + tol)
    BooleanShape{T,OP,TreeOperand{T},TreeOperand{T}}(left, right, place, left_aabb, right_aabb)
end

#---Constructor------------------------------------------------------------------------------------
function BooleanUnion(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanUnion{T}(left, right, place)
end
function BooleanSubtraction(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanSubtraction{T}(left, right, place)
end
function BooleanIntersection(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanIntersection{T}(left, right, place)
end

function BooleanSubtraction(left::NoShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    return left
end

function BooleanIntersection(left::NoShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    return left
end

function BooleanIntersection(left::AbstractShape{T}, right::NoShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    return right
end

#---Utilities---------------------------------------------------------------------------------------

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, shape::BooleanUnion{T}) where {T}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanUnion",(left=left, right=right, placement=placement))
end
function Base.show(io::IO, shape::BooleanSubtraction{T}) where {T}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanSubtraction",(left=left, right=right, placement=placement))
end
function Base.show(io::IO, shape::BooleanIntersection{T}) where {T}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanIntersection",(left=left, right=right, placement=placement))
end


function GeometryBasics.mesh(shape::BooleanUnion{T}; facets::Integer=64) where {T}
    (; left, right, transformation) = shape
    rmesh = mesh(right; facets)
    merge([mesh(left; facets), Mesh(map(c -> _gbpoint(Point3{T}(c) * transformation), coordinates(rmesh)), faces(rmesh))])
end
function GeometryBasics.mesh(shape::BooleanSubtraction{T}; facets::Integer=64) where {T}
    (; left, right, transformation) = shape
    rmesh = mesh(right; facets)
    merge([mesh(left; facets), Mesh(map(c -> _gbpoint(Point3{T}(c) * transformation), coordinates(rmesh)), faces(rmesh))])
end
function GeometryBasics.mesh(shape::BooleanIntersection{T}; facets::Integer=64) where {T}
    (; left, right, transformation) = shape
    rmesh = mesh(right; facets)
    merge([mesh(left; facets), Mesh(map(c -> _gbpoint(Point3{T}(c) * transformation), coordinates(rmesh)), faces(rmesh))])
end

#---Basic functions---------------------------------------------------------------------------------
function extent(shape::BooleanUnion{T})::Tuple{Point3{T},Point3{T}} where {T}
    (; left, right, transformation) = shape
    minLeft, maxLeft = extent(left)
    minRight, maxRight = transform_extent(extent(right),transformation)
    (min.(minLeft, minRight), max.(maxLeft, maxRight))
end
function extent(shape::BooleanSubtraction{T})::Tuple{Point3{T},Point3{T}} where {T}
    extent(shape.left)
end
function extent(shape::BooleanIntersection{T})::Tuple{Point3{T},Point3{T}} where {T}
     (; left, right, transformation) = shape
    minLeft, maxLeft = extent(left)
    minRight, maxRight = transform_extent(extent(right),transformation)
    (max.(minLeft, minRight), min.(maxLeft, maxRight))
end

 
#---Probes for the points on the surface of both operands--------------------------------------------
#   There the classification depends on the neighbourhood of the point: the operands may only touch
#   (a seam of a union, or coinciding faces of a subtraction) or form an edge of the boolean shape. The
#   normals of the operands give the candidate directions, and the points at a small distance along
#   them are classified. The normals alone are not enough: they are averaged over the faces of edges.
_probe_point(point::Point3{T}, dir::Vector3{T}) where T = point + 2 * tolerance(T) * normalize(dir)

# Is the point at a small distance from `point` along `dir` outside both operands of the union?
function _union_probe_outside(shape::BooleanUnion{T}, point::Point3{T}, dir::Vector3{T}) where T
    dir ⋅ dir > 0 || return false
    probe = _probe_point(point, dir)
    inside(shape.left, probe) == kOutside && inside(shape.right, shape.transformation * probe) == kOutside
end

# The 26 directions to the neighbours of a cell in a cubic grid, and the sum of those along which the
# points next to `point` are outside the union: for a point on its surface, the sum is along the
# normal of a face, and along the bisector of an edge. It is zero for a point inside.
const _NEIGHBOUR_DIRECTIONS = Tuple((i, j, k) for i in -1:1 for j in -1:1 for k in -1:1 if (i, j, k) != (0, 0, 0))
function _union_outside_directions(shape::BooleanUnion{T}, point::Point3{T}) where T
    sum = zero(Vector3{T})
    for (i, j, k) in _NEIGHBOUR_DIRECTIONS
        dir = Vector3{T}(i, j, k)
        _union_probe_outside(shape, point, dir) && (sum += normalize(dir))
    end
    return sum
end

# Normal of `shape` at a point of its surface, from a plane through four points of the surface around
# it, found by bisection along the approximate normal `c` (the points are symmetric, so the result is
# exact on spheres and cylinders). Returns zero if the surface is not found.
function _fitted_normal(shape::AbstractShape{T}, point::Point3{T}, c::Vector3{T})::Vector3{T} where T
    δ = min(1000 * tolerance(T), T(0.1))     # the surface is known within the tolerance
    e = abs(c[1]) < T(0.6) ? Vector3{T}(1, 0, 0) : Vector3{T}(0, 1, 0)
    u = normalize(c × e)
    v = c × u
    function surface_point(q)
        lo, hi = -2δ, 2δ                                     # inside at lo, outside at hi
        inside(shape, q + lo * c) == kOutside && return nothing
        inside(shape, q + hi * c) == kOutside || return nothing
        while hi - lo > tolerance(T) / 4
            mid = (lo + hi) / 2
            inside(shape, q + mid * c) == kOutside ? (hi = mid) : (lo = mid)
        end
        return q + (lo + hi) / 2 * c
    end
    pu1, pu2 = surface_point(point + δ * u), surface_point(point - δ * u)
    pv1, pv2 = surface_point(point + δ * v), surface_point(point - δ * v)
    any(isnothing, (pu1, pu2, pv1, pv2)) && return zero(Vector3{T})
    n = (pu1 - pu2) × (pv1 - pv2)
    n ⋅ n > 0 || return zero(Vector3{T})
    n = normalize(n)
    return n ⋅ c < 0 ? -n : n
end

# Normal on an edge of a boolean shape from the approximate normal `c`, refined with the surface around
# the point. One of the `candidates` (normals of the operands) is returned if it agrees, as it is exact.
function _refined_normal(shape::BooleanShape{T}, point::Point3{T}, c::Vector3{T}, candidates::Vector3{T}...)::Vector3{T} where T
    c = normalize(c)
    f = _fitted_normal(shape, point, c)
    f ⋅ f > 0 && (c = f)
    for n in candidates
        n ⋅ n > 0 || continue
        n = normalize(n)
        n ⋅ c > 1 - T(1e-4) && return n
    end
    return c
end

# Is the point at a small distance from `point` along `dir` inside A and outside B?
function _subtraction_probe_inside(shape::BooleanSubtraction{T}, point::Point3{T}, dir::Vector3{T}) where T
    dir ⋅ dir > 0 || return false
    probe = _probe_point(point, dir)
    inside(shape.left, probe) == kInside && inside(shape.right, shape.transformation * probe) == kOutside
end

function inside(shape::BooleanUnion{T}, point::Point3{T})::Int64 where {T}

    course_inside_left = inside(shape.left_aabb, point)
    course_inside_right = inside(shape.right_aabb, point)

    # If the point isn't in either of the BB it is not inside
    !course_inside_left && !course_inside_right && return kOutside
    
    (; left, right, transformation) = shape
    
    # If inside one BB but not the other
    course_inside_left && !course_inside_right && return inside(left, point)
    lpoint = transformation * point
    course_inside_right && !course_inside_left && return inside(right, lpoint)

    positionA = inside(left, point)
    positionA == kInside && return kInside

    positionB = inside(right, lpoint)
    positionB == kInside && return kInside

    if positionA == kSurface && positionB == kSurface
        # On the surface of both operands: on the surface of the union if there are points outside
        # both operands next to it, and inside if the operands touch there (a seam -)(-). The normals of
        # the operands are averaged at their edges, so their sign alone does not tell: the candidate
        # directions are probed instead.
        normalA = normal(left, point)
        normalB = normal(right, lpoint) * transformation
        _union_probe_outside(shape, point, normalA + normalB) && return kSurface
        _union_probe_outside(shape, point, normalA) && return kSurface
        _union_probe_outside(shape, point, normalB) && return kSurface
        return _union_outside_directions(shape, point) == zero(Vector3{T}) ? kInside : kSurface
    elseif positionA == kSurface || positionB == kSurface
        return kSurface
    else
        return kOutside
    end
end

function inside(shape::BooleanIntersection{T}, point::Point3{T})::Int64  where {T}
    
    course_inside_left = inside(shape.left_aabb, point)
    course_inside_right = inside(shape.right_aabb, point)

    # If the point isn't in both BB it is not in the shape
    !(course_inside_left && course_inside_right) && return kOutside

    (; left, right, transformation) = shape
    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kOutside && return kOutside

    positionB = inside(right, lpoint)
    if positionA == kInside && positionB == kInside
        return kInside
    elseif  positionA == kInside && positionB == kSurface ||
            positionA == kSurface && positionB == kInside ||
            positionA == kSurface && positionB == kSurface
        return kSurface
    else
        return kOutside
    end 
end

function inside(shape::BooleanSubtraction{T}, point::Point3{T})::Int64  where {T}
    course_inside_left = inside(shape.left_aabb,point)
    course_inside_right = inside(shape.right_aabb,point)
    (; left, right, transformation) = shape
    # If the point is the left BB but not the right 
    course_inside_left && !course_inside_right && return inside(left, point)

    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kOutside && return kOutside

    positionB = inside(right, lpoint)

    if positionA == kInside && positionB == kOutside
        return kInside;
    elseif positionA == kInside && positionB == kSurface ||
           positionB == kOutside && positionA == kSurface
        return kSurface
    elseif positionA == kSurface && positionB == kSurface
        # On the surface of both operands: on an edge of the subtraction if there are points of A - B
        # next to it, and outside if the faces of A and B coincide. As for the union, the candidate
        # directions given by the normals are probed.
        normalA = normal(left, point)
        normalB = normal(right, lpoint) * transformation
        _subtraction_probe_inside(shape, point, normalB - normalA) && return kSurface
        _subtraction_probe_inside(shape, point, -normalA) && return kSurface
        _subtraction_probe_inside(shape, point, normalB) && return kSurface
        return kOutside
    else
        return kOutside
    end
end

#---Normals------------------------------------------------------------------------------------------
#   On the surface of a single operand the normal is the one of the operand. On the surface of both
#   (edges and seams) the combination of their normals is only an estimate, refined by probing the
#   shape around the point (see _refined_normal).
_unit_or_zero(n::Vector3{T}) where T = n ⋅ n > 0 ? normalize(n) : zero(Vector3{T})

function normal(shape::BooleanUnion{T}, point::Point3{T})::Vector3{T} where {T}
    (; left, right, transformation) = shape
    lpoint = transformation * point
    posA, posB = inside(left, point), inside(right, lpoint)
    if posA == kSurface && posB == kOutside
        return normal(left, point)
    elseif posB == kSurface && posA == kOutside
        return normal(right, lpoint) * transformation
    elseif posA == kSurface && posB == kSurface
        nA, nB = normal(left, point), normal(right, lpoint) * transformation
        nA ⋅ nB > 1 - √eps(T) && return _unit_or_zero(nA + nB)     # coinciding faces
        # Seams and edges: the normals of the operands are averaged over their own edges, and their sum
        # is not the normal of the union when several operands meet (e.g. at a face divided into four).
        # The directions to the outside give an estimate, refined with the surface around the point.
        # An operand normal that agrees with it is returned as it is exact.
        c = _union_outside_directions(shape, point)
        c ⋅ c > 0 || return zero(Vector3{T})                       # on a seam inside the union
        return _refined_normal(shape, point, c, nA + nB, nA, nB)
    end
    return zero(Vector3{T})
end

function normal(shape::BooleanSubtraction{T}, point::Point3{T})::Vector3{T} where {T}
    (; left, right, transformation) = shape
    lpoint = transformation * point
    posA, posB = inside(left, point), inside(right, lpoint)
    if posA == kSurface && posB == kOutside
        return normal(left, point)
    elseif posB == kSurface && posA == kInside
        return -(normal(right, lpoint) * transformation)
    elseif posA == kSurface && posB == kSurface      # edge where the cut meets the outer surface
        nA, nB = normal(left, point), normal(right, lpoint) * transformation
        c = nA - nB
        c ⋅ c > 0 || return zero(Vector3{T})
        return _refined_normal(shape, point, c, nA - nB, nA, -nB)
    end
    return zero(Vector3{T})
end

function normal(shape::BooleanIntersection{T}, point::Point3{T})::Vector3{T} where {T}
    (; left, right, transformation) = shape
    lpoint = transformation * point
    posA, posB = inside(left, point), inside(right, lpoint)
    if posA == kSurface && posB == kInside
        return normal(left, point)
    elseif posB == kSurface && posA == kInside
        return normal(right, lpoint) * transformation
    elseif posA == kSurface && posB == kSurface
        nA, nB = normal(left, point), normal(right, lpoint) * transformation
        c = nA + nB
        c ⋅ c > 0 || return zero(Vector3{T})
        return _refined_normal(shape, point, c, nA + nB, nA, nB)
    end
    return zero(Vector3{T})
end

#---Capacity and surface (no closed form: estimated by sampling, as in Geant4)------------------------
capacity(shape::BooleanShape; npoints::Integer=1_000_000) = estimate_capacity(shape; npoints)
surface(shape::BooleanShape; npoints::Integer=1_000_000) = estimate_surface(shape; npoints)

function safety_to_out(shape::BooleanUnion{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_out(shape::BooleanSubtraction{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_out(shape::BooleanIntersection{T}, point::Point3{T})::T where {T}
    return 0
end

function safety_to_in(shape::BooleanUnion{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_in(shape::BooleanSubtraction{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_in(shape::BooleanIntersection{T}, point::Point3{T})::T where {T}
    return 0
end

function distance_to_out(shape::BooleanUnion{T}, point::Point3{T}, dir::Vector3{T})::T where {T}

    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    (; left, right, transformation) = shape
    # If it intersects one BB but not the other just use that one 
    course_intersects_left && !course_intersects_right && return distance_to_out(left, point, dir)
    course_intersects_right && !course_intersects_left && return distance_to_out(right, transformation * point, transformation * dir)

    dist = T(0)
    positionA = inside(left, point)
    if positionA != kOutside  # point inside A
        while(true)
            distA = distance_to_out(left, point, dir)
            dist += (distA > 0 && distA < T(Inf)) ? distA : 0
            dist += push_tolerance(T)
            npoint = point + dist * dir
            lpoint = transformation * npoint
            if inside(right, lpoint) != kOutside # B could be overlapping with A -- and/or connecting A to another part of A
                ldir = transformation * dir
                distB = distance_to_out(right, lpoint, ldir)
                dist += (distB > 0 && distB < T(Inf)) ? distB : 0
                dist += push_tolerance(T)
                npoint = point + dist * dir
                if inside(left, npoint) == kOutside
                    break
                end
            else
                break
            end
        end
        return dist - push_tolerance(T)
    else
        lpoint = transformation * point
        positionB = inside(right, lpoint)
        if positionB != kOutside  # point inside B
            ldir = transformation * dir
            while(true)
                distB = distance_to_out(right, lpoint, ldir)
                dist += (distB > 0 && distB < T(Inf)) ? distB : 0
                dist += push_tolerance(T) # Give a push
                npoint = point + dist * dir
                if inside(left, npoint) != kOutside # A could be overlapping with B -- and/or connecting B to another part of B
                    ldir = transformation * dir
                    distA = distance_to_out(left, npoint, dir)
                    dist += (distA > 0 && distA < T(Inf)) ? distA : 0
                    dist += push_tolerance(T)
                    npoint = point + dist * dir
                    lpoint = transformation * npoint
                    if inside(right, lpoint) == kOutside
                        break
                    end
                else
                    break
                end
            end
            return dist - push_tolerance(T)
        else
            return T(-1)
        end
    end 
end

function distance_to_out(shape::BooleanIntersection{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    (; left, right, transformation) = shape
    distA = distance_to_out(left, point, dir)
    distB = distance_to_out(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distance_to_out(shape::BooleanSubtraction{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    (; left, right, transformation) = shape
    distA = distance_to_out(left, point, dir)
    distB = distance_to_in(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distance_to_in(shape::BooleanUnion{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    # If it doesn't intersect either BB skip everything
    !course_intersects_left && !course_intersects_right && return T(Inf)
    
    (; left, right, transformation) = shape
    
    # If it intersects one BB but not the other just use that one 
    course_intersects_left && !course_intersects_right && return distance_to_in(left, point, dir)
    course_intersects_right && !course_intersects_left && return distance_to_in(right, transformation * point, transformation * dir)
    
    # If it gets here we have to calculate the distance the complicated way
    distA = distance_to_in(left, point, dir)
    distB = distance_to_in(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distance_to_in(shape::BooleanIntersection{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    # If it doesn't intersect both BB skip everything
    !(course_intersects_left && course_intersects_right) && return T(Inf)

    # If it gets here we have to calculate the distance the complicated way
    (; left, right, transformation) = shape

    positionA = inside(left, point)
    lpoint = transformation * point
    positionB = inside(right, lpoint)

    inLeft = positionA == kInside
    inRight = positionB == kInside

    inLeft && inRight && return T(-1)

    dist = T(0)
    npoint = point
    ldir = transformation * dir
    #  main loop
    while true
        d1 = d2 = zero(T)
        if !inLeft
            d1 = distance_to_in(left, npoint, dir)
            d1 = max(d1, tolerance(T))
            d1 == T(Inf) && return T(Inf)
        end
        if !inRight
            d2 = distance_to_in(right, lpoint, ldir)
            d2 = max(d2, tolerance(T))
            d2 == T(Inf) && return T(Inf)
        end
        if d1 > d2
            # propagate to left shape
            dist += d1
            inLeft = true
            npoint += d1 * dir
            lpoint = transformation * npoint
            # check if propagated point is inside right shape, check is done with a little push
            inRight = inside(right, lpoint + tolerance(T) * ldir) == kInside
            inRight && return dist
            # here inleft=true, inright=false
        else
            # propagate to right shape
            dist += d2
            inRight = true
            # check if propagated point is inside left shape, check is done with a little push
            npoint += d2 * dir
            lpoint = transformation * npoint
            inLeft = inside(left, npoint + tolerance(T) * dir) == kInside
            inLeft && return dist
        end
    end
    return dist
end

function distance_to_in(shape::BooleanSubtraction{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    

    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    # If it doesn't intersect the left shape BB we can skip everyingthing
    !course_intersects_left && return T(Inf)

    (; left, right, transformation) = shape

    # If it doesn't intersect the right we can us just the left distance function 
    !course_intersects_right && return distance_to_in(left, point, dir)

    # If it gets here we have to calculate the distance the complicated way
    lpoint = transformation * point
    ldir = transformation * dir
    positionB = inside(right, lpoint)
    inRight = positionB == kInside

    # `offset` is the distance from `point` to the probe point `npoint`; distances computed at the
    # probe point are measured from there. Probe points are pushed beyond the boundaries crossed.
    npoint = point
    offset = T(0)

    while true
        if inRight
            # propagate to outside of '- / RightShape'
            d1 = distance_to_out(right, lpoint, ldir)
            exit = offset + ((d1 >= 0 && d1 < T(Inf)) ? d1 : zero(T))
            offset = exit + push_tolerance(T)
            npoint = point + offset * dir
            lpoint = transformation * npoint
            # now master outside 'B'; check if inside 'A'
            inside(left, npoint) == kInside && distance_to_out(left, npoint, dir) > push_tolerance(T) && return exit
        end

        # if outside of both we do a max operation master outside '-' and outside '+' ;  find distances to both
        d2 = distance_to_in(left, npoint, dir)
        d2 = max(d2, 0)
        d2 == T(Inf) && return T(Inf)

        d1 = distance_to_in(right, lpoint, ldir)
        if d2 < d1 - tolerance(T)
          return offset + d2
        end

        #   propagate to '-'
        offset += (d1 >= 0 && d1 < T(Inf)) ? d1 + push_tolerance(T) : push_tolerance(T)
        npoint = point + offset * dir
        lpoint = transformation * npoint
        inRight = true
    end
end

