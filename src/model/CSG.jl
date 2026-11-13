#---Boolean shapes of a flat GeometryModel as CSG expressions---------------------------------------------
#   A boolean shape (unions, subtractions, intersections and multi-unions, nested in any way) is stored
#   as a list of primitive components and a postfix program. The algorithms are generic and do not
#   recurse:
#   - inside: evaluate the program with a stack of inside/surface/outside values (2 bits per entry)
#   - distances: step from one component boundary to the next along the ray, and evaluate the program
#     at a point slightly beyond each boundary
#   - normal: normal of a component on whose surface the point is, checked by probing both sides
#   Each primitive function is called from a few places only, so the size of the generated code does
#   not depend on the nesting of the boolean shapes.

"A boolean shape of a GeometryModel (view of a `FlatCSG` record)."
struct CSGShape{T,M} <: AbstractShape{T}
    model::M
    csg::FlatCSG{T}
end

const CSG_MAX_STACK = 32       # 2 bits per entry in a UInt64

@inline _component(c::CSGShape, k) = @inbounds c.model.csg_components[c.csg.first_comp + k + 1]   # k is 0-based
@inline _primitive(f, c::CSGShape, comp::CSGComponent, args...) = dispatch(f, c.model, comp.ref, Val(0), args...)

@inline function _combine(op::Int32, a::UInt64, b::UInt64)::UInt64
    I, S, O = UInt64(kInside), UInt64(kSurface), UInt64(kOutside)
    if op == CSG_UNION
        (a == I || b == I) ? I : (a == S || b == S) ? S : O
    elseif op == CSG_INTERSECTION
        (a == O || b == O) ? O : (a == I && b == I) ? I : S
    else  # CSG_SUBTRACTION
        (a == O || b == I) ? O : (a == I && b == O) ? I : S
    end
end

#---Fast path for subtractions and intersections of two components (A - B, A ∩ B)--------------------
#   From a point inside the shape, the first boundary crossed is the exit: the distance to out is
#   min(distance_to_out(A), distance_to_in(B)) for a subtraction and min(distance_to_out(A),
#   distance_to_out(B)) for an intersection, as in Geant4 and VecGeom, with no probe points. The loop
#   (over the run-time number of components, so not unrolled) keeps a single call site of each
#   primitive function: the GPU compilers inline all of them, and Metal takes much longer to compile
#   a kernel when they are repeated.
@inline _pair_op(c::CSGShape) = @inbounds c.model.csg_program[c.csg.first_op + 3]
@inline _has_pair_exit(c::CSGShape) = c.csg.ncomp == 2 && _pair_op(c) != CSG_UNION

# The point is inside the shape (otherwise the result is -1 or meaningless, as for the primitives)
function _pair_distance_to_out(c::CSGShape{T}, point::Point3{T}, dir::Vector3{T})::T where T
    subtraction = _pair_op(c) == CSG_SUBTRACTION
    invdir = inv.(dir)
    dist = T(Inf)
    for k in 0:c.csg.ncomp-1
        comp = _component(c, k)
        if k == 1 && subtraction         # distance to enter B, clamped to 0 on its surface moving in
            (inside(comp.aabb, point) || intersect(comp.aabb, point, dir, invdir)) || continue
            d = max(_primitive(distance_to_in, c, comp, comp.transformation * point, comp.transformation * dir), zero(T))
        else
            d = _primitive(distance_to_out, c, comp, comp.transformation * point, comp.transformation * dir)
            d < 0 && return T(-1)
        end
        dist = min(dist, d)
    end
    return dist
end

function inside(c::CSGShape{T}, point::Point3{T})::Int64 where T
    inside(c.csg.aabb, point) || return kOutside
    stack = UInt64(0)
    for i in 1:c.csg.nops
        op = @inbounds c.model.csg_program[c.csg.first_op + i]
        if op > 0
            comp = _component(c, op - 1)
            v = inside(comp.aabb, point) ? _primitive(inside, c, comp, comp.transformation * point) : kOutside
            stack = (stack << 2) | UInt64(v)
        else
            b = stack & 0x3; stack >>= 2
            a = stack & 0x3; stack >>= 2
            stack = (stack << 2) | _combine(op, a, b)
        end
    end
    return Int64(stack & 0x3)
end

# Distance from `point` to the nearest boundary of a component ahead along `dir` (Inf if none)
function _next_boundary(c::CSGShape{T}, point::Point3{T}, dir::Vector3{T}, invdir::Vector3{T}) where T
    tmin = T(Inf)
    for k in 0:c.csg.ncomp-1
        comp = _component(c, k)
        (inside(comp.aabb, point) || intersect(comp.aabb, point, dir, invdir)) || continue
        lp, ld = comp.transformation * point, comp.transformation * dir
        pos = _primitive(inside, c, comp, lp)
        entering = pos == kOutside || (pos == kSurface && _primitive(normal, c, comp, lp) ⋅ ld < 0)
        d = pos == kInside ? _primitive(distance_to_out, c, comp, lp, ld) :
            entering && pos == kSurface ? _primitive(distance_to_out, c, comp, lp, ld) :
                                          _primitive(distance_to_in, c, comp, lp, ld)
        d > 0 && d < tmin && (tmin = d)
    end
    return tmin
end

# Classify a probe point beyond a boundary; with a grazing crossing it can still be within the
# tolerance of a surface (at shallow angles the band along the ray is wide, mostly in single
# precision): push it further along the ray until the classification is decisive.
@inline function _classify_probe(c::CSGShape{T}, point, dir, offset) where T
    probe = point + offset * dir
    pos = inside(c, probe)
    push = push_tolerance(T)
    for _ in 1:6                         # pushes doubling up to 32× the first: about 63 in total
        pos == kSurface || break
        offset += push
        push += push
        probe = point + offset * dir
        pos = inside(c, probe)
    end
    return pos, offset
end

function distance_to_in(c::CSGShape{T}, point::Point3{T}, dir::Vector3{T})::T where T
    invdir = inv.(dir)
    (inside(c.csg.aabb, point) || intersect(c.csg.aabb, point, dir, invdir)) || return T(Inf)
    pos = inside(c, point)
    pos == kInside && return T(-1)
    push = push_tolerance(T)
    # on the surface: entering if the point moved along `dir` gets inside (pushed until decisive,
    # a single push can stay within the surface tolerance when leaving at an angle)
    pos == kSurface && first(_classify_probe(c, point, dir, push)) != kOutside && return zero(T)
    offset = zero(T)                     # distance to the probe point (pushed beyond the boundaries)
    probe = point
    for _ in 1:4 * c.csg.ncomp + 16
        t = _next_boundary(c, probe, dir, invdir)
        t == T(Inf) && return T(Inf)
        boundary = offset + t
        pos, offset = _classify_probe(c, point, dir, boundary + push)
        pos != kOutside && return boundary
        probe = point + offset * dir
    end
    return T(Inf)
end

function distance_to_out(c::CSGShape{T}, point::Point3{T}, dir::Vector3{T})::T where T
    _has_pair_exit(c) && return _pair_distance_to_out(c, point, dir)
    inside(c, point) == kOutside && return T(-1)
    invdir = inv.(dir)
    push = push_tolerance(T)
    offset = zero(T)
    probe = point
    boundary = zero(T)
    for _ in 1:4 * c.csg.ncomp + 16
        t = _next_boundary(c, probe, dir, invdir)
        t == T(Inf) && return boundary
        boundary = offset + t
        pos, offset = _classify_probe(c, point, dir, boundary + push)
        pos == kOutside && return boundary
        probe = point + offset * dir
    end
    return boundary
end

function normal(c::CSGShape{T}, point::Point3{T})::Vector3{T} where T
    ε = 2 * tolerance(T)
    for k in 0:c.csg.ncomp-1
        comp = _component(c, k)
        inside(comp.aabb, point) || continue
        lp = comp.transformation * point
        _primitive(inside, c, comp, lp) == kSurface || continue
        n = _primitive(normal, c, comp, lp) * comp.transformation
        n ⋅ n > 0 || continue
        out, in = inside(c, point + ε * n), inside(c, point - ε * n)
        out == kOutside && in != kOutside && return n
        in == kOutside && out != kOutside && return -n      # e.g. the surface of a subtracted shape
    end
    return zero(Vector3{T})
end

extent(c::CSGShape{T}) where T = (c.csg.aabb.min, c.csg.aabb.max)
safety_to_in(c::CSGShape{T}, point::Point3{T}) where T = zero(T)       # conservative
safety_to_out(c::CSGShape{T}, point::Point3{T}) where T = zero(T)
