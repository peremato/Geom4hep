#---Polyhedra---------------------------------------------------------------------------------------
#   Polygonal analogue of the Polycone, as G4Polyhedra: a sequence of sections along z whose cross
#   sections are regular polygons (an outer one, and optionally an inner one making a hole) with
#   `nsides` sides covering the ϕ range [ϕ₀, ϕ₀ + Δϕ]. As in Geant4 and GDML, the radii given at each
#   z plane are the distances from the z axis to the sides (apothems), not to the corners. The sides
#   of a section are planes (tilted if the radius changes along z). The corners are at the angles
#   ϕ₀ + kα (α = Δϕ / nsides), so side k (0-based) has its outward normal at the angle ϕ₀ + (k + ½)α.
#
#   Within the ϕ range a point is in the sector of one side, and the projection of the point on the
#   normal of that side is the largest of all sides: a section is the set of points with
#   z₁ ≤ z ≤ z₂, inside the ϕ range, with projections ≤ rmax(z) on all sides (convex) and not with
#   projections < rmin(z) on all sides (the hole, convex). The intersection of a line with a section
#   is then an interval (convex part) minus at most two intervals: the hole, and the missing ϕ wedge
#   when Δϕ > π (the wedge is convex when Δϕ ≤ π and is clipped together with the sides).
#---------------------------------------------------------------------------------------------------
"Section of a `Polyhedra` between `z1` and `z2` (radii are distances to the sides)."
struct PolyhedraSection{T<:AbstractFloat}
    z1::T
    z2::T
    rmin1::T
    rmax1::T
    rmin2::T
    rmax2::T
    smin::T          # slopes of the radii, d(radius)/dz
    smax::T
    kmin::T          # 1/√(1 + slope²): normalization of the equations of the side planes
    kmax::T
end
function PolyhedraSection{T}(z1, z2, rmin1, rmax1, rmin2, rmax2) where T<:AbstractFloat
    h = z2 - z1
    smin, smax = (rmin2 - rmin1) / h, (rmax2 - rmax1) / h
    PolyhedraSection{T}(z1, z2, rmin1, rmax1, rmin2, rmax2, smin, smax, 1 / sqrt(1 + smin^2), 1 / sqrt(1 + smax^2))
end
@inline _hollow(s::PolyhedraSection) = s.rmin1 > 0 || s.rmin2 > 0

struct Polyhedra{T<:AbstractFloat, VS<:AbstractVector{PolyhedraSection{T}}} <: AbstractShape{T}
    sections::VS       # sections (non-zero length) in increasing z
    nsides::Int32      # number of sides in the ϕ range
    ϕ₀::T              # starting ϕ (radians)
    Δϕ::T              # ϕ extent (radians)
    fullϕ::Bool
    wedge::Wedge{T}
end

"""
    Polyhedra{T}(nsides, rmin, rmax, z, ϕ₀=0, Δϕ=2π)

Polyhedra (as G4Polyhedra) with `nsides` sides in the ϕ range [ϕ₀, ϕ₀+Δϕ] and z planes `z` with
inner and outer radii `rmin` and `rmax`. The radii are the distances from the z axis to the sides.
"""
function Polyhedra{T}(nsides::Integer, rmin::AbstractVector, rmax::AbstractVector, z::AbstractVector, ϕ₀=0, Δϕ=2π) where T<:AbstractFloat
    L = length(z)
    length(rmin) == length(rmax) == L || error("Polyhedra: vectors rmin, rmax and z have different sizes")
    L >= 2 || error("Polyhedra: at least two z planes are needed")
    if z[end] < z[1]                              # planes given in decreasing z
        rmin, rmax, z = reverse(rmin), reverse(rmax), reverse(z)
    end
    all(i -> z[i+1] >= z[i], 1:L-1) || error("Polyhedra: the z planes are not ordered")
    all(i -> 0 <= rmin[i] <= rmax[i], 1:L) || error("Polyhedra: invalid radii")
    fullϕ = Δϕ >= 2π - 1e-9
    fullϕ && (Δϕ = 2π)
    nsides >= 1 && Δϕ / nsides < π || error("Polyhedra: invalid number of sides ($nsides for Δϕ=$Δϕ)")
    sections = [PolyhedraSection{T}(z[i], z[i+1], rmin[i], rmax[i], rmin[i+1], rmax[i+1]) for i in 1:L-1 if z[i+1] > z[i]]
    N = length(sections)
    Polyhedra{T, SVector{N,PolyhedraSection{T}}}(SVector{N}(sections), nsides, ϕ₀, Δϕ, fullϕ, Wedge{T}(ϕ₀, Δϕ))
end

function Base.show(io::IO, ph::Polyhedra{T}) where T
    z = [[s.z1 for s in ph.sections]; ph.sections[end].z2]
    print(io, "Polyhedra{$T}", (nsides=Int(ph.nsides), nsections=length(ph.sections), z=z, ϕ₀=ph.ϕ₀, Δϕ=ph.Δϕ))
end

nsections(ph::Polyhedra) = length(ph.sections)

#---Sides and constraints----------------------------------------------------------------------------
@inline _side_angle(ph::Polyhedra{T}, k) where T = ph.ϕ₀ + (k + T(0.5)) * ph.Δϕ / ph.nsides

# Side (0-based) whose sector contains the point (x, y); outside the ϕ range, the side at the nearest end
@inline function _side_index(ph::Polyhedra{T}, x, y) where T
    n = ph.nsides
    n == 1 && return Int32(0)
    δ = mod(atan(y, x) - ph.ϕ₀, T(2π))
    if δ < ph.Δϕ
        return min(unsafe_trunc(Int32, δ * n / ph.Δϕ), n - Int32(1))
    else
        return δ - ph.Δϕ < T(2π) - δ ? n - Int32(1) : Int32(0)
    end
end

# Normal (cos, sin) of the side of the sector of (x, y) and the projection of the point on it
@inline function _side_projection(ph::Polyhedra{T}, x, y) where T
    s, c = sincos(_side_angle(ph, _side_index(ph, x, y)))
    return c, s, c * x + s * y
end

# Distance-like value of the ϕ wedge (positive outside the ϕ range)
@inline function _wedge_value(ph::Polyhedra{T}, x, y) where T
    ph.fullϕ && return T(-Inf)
    w = ph.wedge
    d1 = w.normal1[1] * x + w.normal1[2] * y          # inward normals of the ϕ planes
    d2 = w.normal2[1] * x + w.normal2[2] * y
    return ph.Δϕ <= π ? -min(d1, d2) : -max(d1, d2)
end

# Largest of the (normalized) constraints of a section except the z planes: ≤ 0 inside, and a lower
# bound of the distance to the section ignoring its z planes. `proj` is the projection of the point on
# the normal of the side of its sector.
@inline function _radial_value(ph::Polyhedra{T}, sec::PolyhedraSection{T}, proj, x, y, z) where T
    t = z - sec.z1
    v = (proj - (sec.rmax1 + sec.smax * t)) * sec.kmax
    _hollow(sec) && (v = max(v, (sec.rmin1 + sec.smin * t - proj) * sec.kmin))
    return max(v, _wedge_value(ph, x, y))
end
@inline _section_value(ph::Polyhedra, sec, proj, x, y, z) = max(_radial_value(ph, sec, proj, x, y, z), sec.z1 - z, z - sec.z2)

# Index of the section containing z (the first one if on a shared plane), 0 below, N+1 above
@inline function _section_at(ph::Polyhedra, z)
    secs = ph.sections
    N = length(secs)
    z < secs[1].z1 && return 0
    for i in 1:N
        z <= secs[i].z2 && return i
    end
    return N + 1
end

#---Basic functions--------------------------------------------------------------------------------
function capacity(ph::Polyhedra{T})::T where T
    f = ph.nsides * tan(ph.Δϕ / ph.nsides / 2)          # area of the polygon of apothem a: f a²
    sum(ph.sections) do s
        h = s.z2 - s.z1
        f * h / 3 * ((s.rmax1^2 + s.rmax1 * s.rmax2 + s.rmax2^2) - (s.rmin1^2 + s.rmin1 * s.rmin2 + s.rmin2^2))
    end
end

function surface(ph::Polyhedra{T})::T where T
    n = ph.nsides
    α = ph.Δϕ / n
    f = n * tan(α / 2)
    area(a, b) = f * (b^2 - a^2)                           # polygonal annulus
    secs = ph.sections
    N = length(secs)
    surf = zero(T)
    for s in secs
        h = s.z2 - s.z1
        # sides: n trapezoids with parallel sides 2 a tan(α/2) and height √(h² + Δa²)
        surf += f * ((s.rmax1 + s.rmax2) * sqrt(h^2 + (s.rmax2 - s.rmax1)^2) + (s.rmin1 + s.rmin2) * sqrt(h^2 + (s.rmin2 - s.rmin1)^2))
        # ϕ faces: trapezoids between the corners (radius / cos(α/2)) at the two ends
        ph.fullϕ || (surf += h * ((s.rmax1 - s.rmin1) + (s.rmax2 - s.rmin2)) / cos(α / 2))
    end
    # faces perpendicular to z: end caps and the steps between consecutive sections
    surf += area(secs[1].rmin1, secs[1].rmax1) + area(secs[N].rmin2, secs[N].rmax2)
    for i in 1:N-1
        lo, hi = secs[i], secs[i+1]
        a = max(lo.rmin2, hi.rmin1)
        overlap = area(a, max(min(lo.rmax2, hi.rmax1), a))
        surf += area(lo.rmin2, lo.rmax2) + area(hi.rmin1, hi.rmax1) - 2overlap
    end
    return surf
end

function extent(ph::Polyhedra{T})::Tuple{Point3{T},Point3{T}} where T
    # bounding box of the corners of the polygons at all z planes
    n = ph.nsides
    α = ph.Δϕ / n
    xmin = ymin = T(Inf)
    xmax = ymax = T(-Inf)
    for s in ph.sections, a in (s.rmin1, s.rmax1, s.rmin2, s.rmax2), k in 0:n
        sn, c = sincos(ph.ϕ₀ + k * α)
        R = a / cos(α / 2)
        xmin, xmax = min(xmin, R * c), max(xmax, R * c)
        ymin, ymax = min(ymin, R * sn), max(ymax, R * sn)
    end
    return (Point3{T}(xmin, ymin, ph.sections[1].z1), Point3{T}(xmax, ymax, ph.sections[end].z2))
end

function inside(ph::Polyhedra{T}, point::Point3{T})::Int64 where T
    x, y, z = point
    tol = tolerance(T) / 2
    secs = ph.sections
    N = length(secs)
    (z < secs[1].z1 - tol || z > secs[N].z2 + tol) && return kOutside
    _, _, proj = _side_projection(ph, x, y)
    lo, hi = _section_at(ph, z - tol), _section_at(ph, z + tol)
    lo < 1 && (lo = hi)
    hi > N && (hi = lo)
    if lo == hi
        v = _section_value(ph, secs[lo], proj, x, y, z)
        return v < -tol ? kInside : v <= tol ? kSurface : kOutside
    end
    # close to a plane shared by two sections: it is not a boundary where both sections are present
    vlo = _radial_value(ph, secs[lo], proj, x, y, z)
    vhi = _radial_value(ph, secs[hi], proj, x, y, z)
    vlo < -tol && vhi < -tol && return kInside
    vlo > tol && vhi > tol && return kOutside
    return kSurface
end
Base.contains(ph::Polyhedra{T}, p::Point3{T}) where T = inside(ph, p) == kInside

function safety_to_in(ph::Polyhedra{T}, point::Point3{T})::T where T
    x, y, z = point
    _, _, proj = _side_projection(ph, x, y)
    safety = T(Inf)
    for sec in ph.sections
        safety = min(safety, _section_value(ph, sec, proj, x, y, z))
    end
    return max(safety, zero(T))
end

function safety_to_out(ph::Polyhedra{T}, point::Point3{T})::T where T
    x, y, z = point
    pos = inside(ph, point)
    pos == kOutside && return T(-1)
    pos == kSurface && return zero(T)
    _, _, proj = _side_projection(ph, x, y)
    i = clamp(_section_at(ph, z), 1, length(ph.sections))
    return max(-_section_value(ph, ph.sections[i], proj, x, y, z), zero(T))   # conservative at shared planes
end

function normal(ph::Polyhedra{T}, point::Point3{T})::Vector3{T} where T
    inside(ph, point) == kSurface || return zero(Vector3{T})
    x, y, z = point
    tol = tolerance(T) / 2
    secs = ph.sections
    N = length(secs)
    c, s, proj = _side_projection(ph, x, y)
    w = ph.wedge
    d1 = w.normal1[1] * x + w.normal1[2] * y
    d2 = w.normal2[1] * x + w.normal2[2] * y
    n = zero(Vector3{T})
    for i in 1:N
        sec = secs[i]
        sec.z1 - tol <= z <= sec.z2 + tol || continue
        _section_value(ph, sec, proj, x, y, z) <= tol || continue
        t = z - sec.z1
        abs(proj - (sec.rmax1 + sec.smax * t)) * sec.kmax <= tol && (n += Vector3{T}(c, s, -sec.smax) * sec.kmax)
        _hollow(sec) && abs(sec.rmin1 + sec.smin * t - proj) * sec.kmin <= tol && (n += Vector3{T}(-c, -s, sec.smin) * sec.kmin)
        if !ph.fullϕ           # ϕ planes (for Δϕ > π a plane is a boundary only where the other is not inside)
            abs(d1) <= tol && (ph.Δϕ <= π || d2 <= tol) && (n -= Vector3{T}(w.normal1[1], w.normal1[2], 0))
            abs(d2) <= tol && (ph.Δϕ <= π || d1 <= tol) && (n -= Vector3{T}(w.normal2[1], w.normal2[2], 0))
        end
        # z planes, unless the neighbouring section covers the point there
        if abs(z - sec.z1) <= tol && !(i > 1 && _radial_value(ph, secs[i-1], proj, x, y, z) < -tol)
            n -= Vector3{T}(0, 0, 1)
        end
        if abs(z - sec.z2) <= tol && !(i < N && _radial_value(ph, secs[i+1], proj, x, y, z) < -tol)
            n += Vector3{T}(0, 0, 1)
        end
    end
    n ⋅ n > 0 || return zero(Vector3{T})
    return normalize(n)
end

# For a point on the surface, whether moving along `dir` goes into the shape: the distance-like value
# of the shape (smallest of the values of the sections near the point) decreases along `dir`. Cheaper
# than the normal (code size matters on GPUs, where all functions are inlined).
@inline function _union_value(ph::Polyhedra{T}, p::Point3{T}) where T
    x, y, z = p
    tol = tolerance(T)
    _, _, proj = _side_projection(ph, x, y)
    v = T(Inf)
    for sec in ph.sections
        sec.z1 - tol <= z <= sec.z2 + tol && (v = min(v, _section_value(ph, sec, proj, x, y, z)))
    end
    return v
end
@inline _moving_in(ph::Polyhedra{T}, p::Point3{T}, d::Vector3{T}) where T =
    _union_value(ph, p + tolerance(T) * d) < _union_value(ph, p)

#---Intersections of a line with a section------------------------------------------------------------
# Clip the interval (lo, hi) of the line p + t d with the half-space g0 + t g1 ≤ 0
@inline function _clip(lo::T, hi::T, g0::T, g1::T) where T
    if g1 > 0
        hi = min(hi, -g0 / g1)
    elseif g1 < 0
        lo = max(lo, -g0 / g1)
    elseif g0 > 0
        lo, hi = T(Inf), T(-Inf)
    end
    return lo, hi
end

# The (up to three) intervals of the line p + t d inside a section, in increasing t (empty: u ≥ v)
@inline function _section_pieces(ph::Polyhedra{T}, sec::PolyhedraSection{T}, p::Point3{T}, d::Vector3{T}) where T
    lo, hi = _clip(T(-Inf), T(Inf), sec.z1 - p[3], -d[3])
    lo, hi = _clip(lo, hi, p[3] - sec.z2, d[3])
    hollow = _hollow(sec)
    A = sec.rmax1 + sec.smax * (p[3] - sec.z1)        # radii at the z of p
    B = sec.rmin1 + sec.smin * (p[3] - sec.z1)
    c1, d1 = T(-Inf), T(Inf)                          # hole: inside the inner polygon
    for k in 0:ph.nsides-1
        s, c = sincos(_side_angle(ph, k))
        pr, dr = c * p[1] + s * p[2], c * d[1] + s * d[2]
        lo, hi = _clip(lo, hi, pr - A, dr - sec.smax * d[3])
        hollow && ((c1, d1) = _clip(c1, d1, pr - B, dr - sec.smin * d[3]))
    end
    hollow || (c1 = d1 = T(Inf))
    c2 = d2 = T(Inf)                                  # missing ϕ wedge (Δϕ > π)
    if !ph.fullϕ
        w = ph.wedge
        g1, e1 = w.normal1[1] * p[1] + w.normal1[2] * p[2], w.normal1[1] * d[1] + w.normal1[2] * d[2]
        g2, e2 = w.normal2[1] * p[1] + w.normal2[2] * p[2], w.normal2[1] * d[1] + w.normal2[2] * d[2]
        if ph.Δϕ <= π                                 # convex wedge: clip
            lo, hi = _clip(lo, hi, -g1, -e1)
            lo, hi = _clip(lo, hi, -g2, -e2)
        else
            c2, d2 = _clip(T(-Inf), T(Inf), g1, e1)
            c2, d2 = _clip(c2, d2, g2, e2)
        end
    end
    # [lo, hi] minus the holes (c1, d1) and (c2, d2)
    c1 < d1 || ((c1, d1) = (c2, d2); c2 = d2 = T(Inf))
    c2 < d2 || (c2 = d2 = T(Inf))
    c2 < c1 && ((c1, d1, c2, d2) = (c2, d2, c1, d1))
    if c2 <= d1                                       # overlapping holes
        d1 = max(d1, d2)
        c2 = d2 = T(Inf)
    end
    return ((lo, min(hi, c1)), (max(lo, d1), min(hi, c2)), (max(lo, d2), hi))
end

# Range of sections to visit along the ray, in its z order
@inline function _section_range(ph::Polyhedra, d)
    N = length(ph.sections)
    return d[3] >= 0 ? (1, N, 1) : (N, 1, -1)
end

#---Distances-----------------------------------------------------------------------------------------
function distance_to_in(ph::Polyhedra{T}, point::Point3{T}, dir::Vector3{T})::T where T
    pos = inside(ph, point)
    pos == kInside && return T(-1)
    tol = tolerance(T) / 2
    entering = pos == kSurface && _moving_in(ph, point, dir)
    z, dz = point[3], dir[3]
    best = T(Inf)
    first, last, step = _section_range(ph, dir)
    for i in first:step:last
        sec = ph.sections[i]
        # sections behind the point, and beyond the best entry found
        dz > 0 && (sec.z2 < z - tol && continue; sec.z1 > z + best * dz + tol && break)
        dz < 0 && (sec.z1 > z + tol && continue; sec.z2 < z + best * dz - tol && break)
        dz == 0 && !(sec.z1 - tol <= z <= sec.z2 + tol) && continue
        for (u, v) in _section_pieces(ph, sec, point, dir)
            u < v || continue
            v - max(u, zero(T)) > tol || continue                # behind, or only grazing
            pos == kSurface && !entering && u <= 0 && continue    # the piece being left
            u < best && (best = u)
            break
        end
    end
    return best
end

function distance_to_out(ph::Polyhedra{T}, point::Point3{T}, dir::Vector3{T})::T where T
    pos = inside(ph, point)
    pos == kOutside && return T(-1)
    pos == kSurface && !_moving_in(ph, point, dir) && return zero(T)     # leaving
    tol = tolerance(T) / 2
    z, dz = point[3], dir[3]
    t = zero(T)
    started = false
    first, last, step = _section_range(ph, dir)
    for i in first:step:last
        sec = ph.sections[i]
        dz > 0 && (sec.z2 < z - tol && continue; sec.z1 > z + t * dz + tol && break)
        dz < 0 && (sec.z1 > z + tol && continue; sec.z2 < z + t * dz - tol && break)
        dz == 0 && !(sec.z1 - tol <= z <= sec.z2 + tol) && continue
        for (u, v) in _section_pieces(ph, sec, point, dir)
            (u < v && v > t) || continue
            if started
                u <= t + tol && (t = v)          # continues the piece (shared z plane)
            elseif u <= tol || pos == kSurface   # the piece of the point
                t = v
                started = true
            end
        end
    end
    return t
end

#---Drawing functions-------------------------------------------------------------------------------
#   Each section is drawn as a closed polyhedron: corners k = 0..n of the outer and inner polygons at
#   both ends; outer and inner sides, the faces at the two z planes and the ϕ faces.
function GeometryBasics.coordinates(ph::Polyhedra{T}, facets=36) where T<:AbstractFloat
    n = ph.nsides
    α = ph.Δϕ / n
    points = Point3{T}[]
    for sec in ph.sections, (z, rmin, rmax) in ((sec.z1, sec.rmin1, sec.rmax1), (sec.z2, sec.rmin2, sec.rmax2)), a in (rmax, rmin), k in 0:n
        s, c = sincos(ph.ϕ₀ + k * α)
        R = a / cos(α / 2)
        push!(points, Point3{T}(R * c, R * s, z))
    end
    return points
end

function GeometryBasics.faces(ph::Polyhedra{T}, facets=36) where T<:AbstractFloat
    n = ph.nsides
    quads = QuadFace{Int}[]
    for i in eachindex(ph.sections)
        off = (i - 1) * 4 * (n + 1)
        O1(k) = off + k + 1                  # outer and inner corners at z1 and z2
        I1(k) = off + (n + 1) + k + 1
        O2(k) = off + 2(n + 1) + k + 1
        I2(k) = off + 3(n + 1) + k + 1
        for k in 0:n-1
            push!(quads, QuadFace{Int}(O1(k), O1(k+1), O2(k+1), O2(k)))    # outer side
            push!(quads, QuadFace{Int}(I1(k), I2(k), I2(k+1), I1(k+1)))    # inner side
            push!(quads, QuadFace{Int}(O1(k), I1(k), I1(k+1), O1(k+1)))    # bottom
            push!(quads, QuadFace{Int}(O2(k), O2(k+1), I2(k+1), I2(k)))    # top
        end
        if !ph.fullϕ
            push!(quads, QuadFace{Int}(O1(0), O2(0), I2(0), I1(0)))        # ϕ₀ face
            push!(quads, QuadFace{Int}(O1(n), I1(n), I2(n), O2(n)))        # ϕ₀ + Δϕ face
        end
    end
    return quads
end
