#---Torus-------------------------------------------------------------------------------------------
#   Torus segment, as G4Torus: the points at a distance between `rmin` and `rmax` from the circle of
#   radius `rtor` around the z axis (the central circle), within the ϕ range [ϕ₀, ϕ₀ + Δϕ].
#
#   The intersections with rays are not computed with the usual quartic equation: for the thin tori
#   of real detectors (rtor/rmax ~ 100-500) its coefficients cancel catastrophically (most digits are
#   lost in double precision, and all in single precision). Instead, the distance D(t) from the ray to
#   the central circle, √((ρ - rtor)² + z²), is evaluated directly (accurate in Float32) and the ray is
#   advanced with steps that cannot jump over a crossing of D with a radius: D is 1-Lipschitz, and its
#   second derivative along a line is bounded,
#       D'' = (F''/2 - D'²)/D   with F = D², F'' = 2 - 2 rtor ρ'' ∈ [2 - 2 rtor/ρ, 2],
#   so D'' ≥ -rtor/(ρ D) and D'' ≤ 1/D. With these bounds the steps behave as Newton steps close to a
#   crossing and remain long along grazing rays.
#---------------------------------------------------------------------------------------------------
"""
    Torus{T}(rmin, rmax, rtor, ϕ₀=0, Δϕ=2π)

Torus or torus segment (as G4Torus): the points at a distance between `rmin` and `rmax` from the
circle of radius `rtor` in the xy plane (the swept radius), within the ϕ range [`ϕ₀`, `ϕ₀ + Δϕ`].
"""
struct Torus{T<:AbstractFloat} <: AbstractShape{T}
    rmin::T            # inner radius of the tube (0 for a solid tube)
    rmax::T            # outer radius of the tube
    rtor::T            # radius of the central circle
    ϕ₀::T              # starting ϕ (radians)
    Δϕ::T              # ϕ extent (radians)
    fullϕ::Bool
    wedge::Wedge{T}
    function Torus{T}(rmin, rmax, rtor, ϕ₀=0, Δϕ=2π) where T<:AbstractFloat
        0 <= rmin < rmax || error("Torus: invalid radii rmin=$rmin, rmax=$rmax")
        rmax < rtor || error("Torus: rmax=$rmax must be smaller than rtor=$rtor")
        Δϕ > 0 || error("Torus: invalid Δϕ=$Δϕ")
        fullϕ = Δϕ >= 2π - 1e-9
        fullϕ && (Δϕ = 2π)
        new(rmin, rmax, rtor, ϕ₀, Δϕ, fullϕ, Wedge{T}(ϕ₀, Δϕ))
    end
end

function Base.show(io::IO, tor::Torus{T}) where T
    print(io, "Torus{$T}", (rmin=tor.rmin, rmax=tor.rmax, rtor=tor.rtor, ϕ₀=tor.ϕ₀, Δϕ=tor.Δϕ))
end

#---Basic functions--------------------------------------------------------------------------------
capacity(tor::Torus) = π * (tor.rmax^2 - tor.rmin^2) * tor.rtor * tor.Δϕ           # Pappus
function surface(tor::Torus{T}) where T
    area = 2π * (tor.rmax + tor.rmin) * tor.rtor * tor.Δϕ
    tor.fullϕ || (area += 2π * (tor.rmax^2 - tor.rmin^2))                             # two ϕ faces
    return area
end
function extent(tor::Torus{T}) where T
    lo, hi = extent(Tube{T}(tor.rtor - tor.rmax, tor.rtor + tor.rmax, tor.rmax, tor.ϕ₀, tor.Δϕ))
    return (lo, hi)
end

# Distance from q to the central circle, and its rate along the unit vector d (and ρ of q)
@inline function _circle_distance(tor::Torus{T}, q, d) where T
    ρ = sqrt(q[1]^2 + q[2]^2)
    u = ρ - tor.rtor
    D = sqrt(u^2 + q[3]^2)
    dρ = ρ > 0 ? (q[1] * d[1] + q[2] * d[2]) / ρ : sqrt(d[1]^2 + d[2]^2)
    dD = D > 0 ? (u * dρ + q[3] * d[3]) / D : one(T)
    return ρ, D, dD
end

# Distance-like value of the ϕ wedge (positive outside the ϕ range)
@inline function _wedge_value(tor::Torus{T}, x, y) where T
    tor.fullϕ && return T(-Inf)
    w = tor.wedge
    d1 = w.normal1[1] * x + w.normal1[2] * y          # inward normals of the ϕ planes
    d2 = w.normal2[1] * x + w.normal2[2] * y
    return tor.Δϕ <= π ? -min(d1, d2) : -max(d1, d2)
end

# Largest of the constraints (≤ 0 inside): outer and inner surfaces of the tube and the ϕ range
@inline function _torus_value(tor::Torus{T}, p) where T
    ρ = sqrt(p[1]^2 + p[2]^2)
    D = sqrt((ρ - tor.rtor)^2 + p[3]^2)
    v = D - tor.rmax
    tor.rmin > 0 && (v = max(v, tor.rmin - D))
    return max(v, _wedge_value(tor, p[1], p[2])), D
end

function inside(tor::Torus{T}, point::Point3{T})::Int64 where T
    v, _ = _torus_value(tor, point)
    tol = tolerance(T) / 2
    return v < -tol ? kInside : v <= tol ? kSurface : kOutside
end
Base.contains(tor::Torus{T}, p::Point3{T}) where T = inside(tor, p) == kInside

function safety_to_in(tor::Torus{T}, point::Point3{T})::T where T
    v, _ = _torus_value(tor, point)                   # lower bound of the distance (exact if full ϕ)
    return max(v, zero(T))
end

function safety_to_out(tor::Torus{T}, point::Point3{T})::T where T
    v, _ = _torus_value(tor, point)
    v > tolerance(T) / 2 && return T(-1)
    return max(-v, zero(T))
end

function normal(tor::Torus{T}, point::Point3{T})::Vector3{T} where T
    inside(tor, point) == kSurface || return zero(Vector3{T})
    x, y, z = point
    tol = tolerance(T) / 2
    ρ = sqrt(x^2 + y^2)
    u = ρ - tor.rtor
    D = sqrt(u^2 + z^2)
    n = zero(Vector3{T})
    if D > 0 && ρ > 0
        g = Vector3{T}(u * x / ρ, u * y / ρ, z) / D   # gradient of D (away from the central circle)
        abs(D - tor.rmax) <= tol && (n += g)
        tor.rmin > 0 && abs(D - tor.rmin) <= tol && (n -= g)
    end
    if !tor.fullϕ          # ϕ planes (for Δϕ > π a plane is a boundary only where the other is not inside)
        w = tor.wedge
        d1 = w.normal1[1] * x + w.normal1[2] * y
        d2 = w.normal2[1] * x + w.normal2[2] * y
        abs(d1) <= tol && (tor.Δϕ <= π || d2 <= tol) && (n -= Vector3{T}(w.normal1[1], w.normal1[2], 0))
        abs(d2) <= tol && (tor.Δϕ <= π || d1 <= tol) && (n -= Vector3{T}(w.normal2[1], w.normal2[2], 0))
    end
    n ⋅ n > 0 || return zero(Vector3{T})
    return normalize(n)
end

# For a point on the surface, whether moving along `d` goes into the shape: the constraints that are
# active (the surfaces the point is on) all decrease along `d` (rates computed analytically)
@inline function _moving_in(tor::Torus{T}, p::Point3{T}, d::Vector3{T}) where T
    tol = tolerance(T) / 2
    _, D, dD = _circle_distance(tor, p, d)
    wv, wr = T(-Inf), T(-Inf)
    if !tor.fullϕ
        w = tor.wedge
        g1, e1 = w.normal1[1] * p[1] + w.normal1[2] * p[2], w.normal1[1] * d[1] + w.normal1[2] * d[2]
        g2, e2 = w.normal2[1] * p[1] + w.normal2[2] * p[2], w.normal2[1] * d[1] + w.normal2[2] * d[2]
        wv, wr = (tor.Δϕ <= π) == (g1 <= g2) ? (-g1, -e1) : (-g2, -e2)
    end
    constraints = ((D - tor.rmax, dD), tor.rmin > 0 ? (tor.rmin - D, -dD) : (T(-Inf), T(-Inf)), (wv, wr))
    vmax = maximum(first, constraints)
    rate = T(-Inf)
    for (v, r) in constraints
        v >= vmax - tol && (rate = max(rate, r))
    end
    return rate < 0
end

#---Intersections with rays---------------------------------------------------------------------------
# First t' ≥ t (up to tend) where the distance D to the central circle crosses the radius c, going
# down (D > c before) or up (D < c before). The steps cannot jump over a crossing (see the header).
# Returns Inf if there is none; when searching up, if the iterations are exhausted, the point reached.
@inline function _crossing(tor::Torus{T}, p::Point3{T}, d::Vector3{T}, t::T, tend::T, c::T, down::Bool) where T
    R = tor.rtor
    ε = tolerance(T) / 4
    μ = sqrt(2ε * c)                  # step along a ray tangent to the surface (within ε)
    n = 0
    @nounroll while t <= tend && n < 200
        ρ, D, dD = _circle_distance(tor, p + t * d, d)
        if down
            g = D - c
            g <= ε && dD < 0 && return t
            K = 2R / (ρ * c)                                        # -D'' ≤ K for steps ≤ ρ/2
            s = ρ > 0 ? min((dD + sqrt(dD^2 + 2K * max(g, zero(T)))) / K, ρ / 2) : zero(T)
        else
            g = c - D
            g <= ε && dD > 0 && return t
            K = 2 / D                                               # D'' ≤ K for steps ≤ D/2
            s = D > 0 ? min((-dD + sqrt(dD^2 + 2K * max(g, zero(T)))) / K, D / 2) : zero(T)
        end
        s = max(s, g)                                              # D is 1-Lipschitz
        g <= ε && (s = max(s, μ))
        t += s
        n += 1
    end
    return down || t > tend ? T(Inf) : t
end

# Interval of the ray inside the bounding region of the torus: |z| ≤ rmax, ρ ≤ rtor + rmax (with margins)
@inline function _torus_bounds(tor::Torus{T}, p::Point3{T}, d::Vector3{T}) where T
    r = tor.rmax * T(1.05) + tolerance(T)
    lo, hi = _clip(T(-Inf), T(Inf), -p[3] - r, -d[3])
    lo, hi = _clip(lo, hi, p[3] - r, d[3])
    Rc = tor.rtor + r
    a = d[1]^2 + d[2]^2
    b = p[1] * d[1] + p[2] * d[2]
    c0 = p[1]^2 + p[2]^2 - Rc^2
    if a > 0
        disc = b^2 - a * c0
        disc < 0 && return T(Inf), T(-Inf)
        sq = sqrt(disc)
        lo, hi = max(lo, (-b - sq) / a), min(hi, (-b + sq) / a)
    elseif c0 > 0
        return T(Inf), T(-Inf)
    end
    return lo, hi
end

# The (up to two) intervals of the ray within the ϕ range (empty: a ≥ b)
@inline function _wedge_intervals(tor::Torus{T}, p::Point3{T}, d::Vector3{T}) where T
    tor.fullϕ && return (T(-Inf), T(Inf), T(Inf), T(-Inf))
    w = tor.wedge
    g1, e1 = w.normal1[1] * p[1] + w.normal1[2] * p[2], w.normal1[1] * d[1] + w.normal1[2] * d[2]
    g2, e2 = w.normal2[1] * p[1] + w.normal2[2] * p[2], w.normal2[1] * d[1] + w.normal2[2] * d[2]
    if tor.Δϕ <= π                                   # convex wedge
        lo, hi = _clip(T(-Inf), T(Inf), -g1, -e1)
        lo, hi = _clip(lo, hi, -g2, -e2)
        return (lo, hi, T(Inf), T(-Inf))
    else                                             # outside of the missing (convex) wedge
        c, e = _clip(T(-Inf), T(Inf), g1, e1)
        c, e = _clip(c, e, g2, e2)
        return c < e ? (T(-Inf), c, e, T(Inf)) : (T(-Inf), T(Inf), T(Inf), T(-Inf))
    end
end

function distance_to_out(tor::Torus{T}, point::Point3{T}, dir::Vector3{T})::T where T
    pos = inside(tor, point)
    pos == kOutside && return T(-1)
    pos == kSurface && !_moving_in(tor, point, dir) && return zero(T)
    tol = tolerance(T) / 2
    a1, b1, a2, b2 = _wedge_intervals(tor, point, dir)
    tend = a1 <= tol && b1 > 0 ? b1 : a2 <= tol && b2 > 0 ? b2 : zero(T)      # leaving the ϕ range
    tend = min(tend, _crossing(tor, point, dir, zero(T), tend, tor.rmax, false))
    tor.rmin > 0 && (tend = min(tend, _crossing(tor, point, dir, zero(T), tend, tor.rmin, true)))
    return tend
end

function distance_to_in(tor::Torus{T}, point::Point3{T}, dir::Vector3{T})::T where T
    pos = inside(tor, point)
    pos == kInside && return T(-1)
    pos == kSurface && _moving_in(tor, point, dir) && return zero(T)
    tn, tf = _torus_bounds(tor, point, dir)
    (tn < tf && tf >= 0) || return T(Inf)
    a1, b1, a2, b2 = _wedge_intervals(tor, point, dir)
    tol = tolerance(T) / 2
    ε = tolerance(T) / 4
    t = max(tn, zero(T))
    n = 0
    @nounroll while n < 8
        n += 1
        # where the ray is with respect to the tube (at its surfaces, according to the direction)
        _, D, dD = _circle_distance(tor, point + t * dir, dir)
        outside = D > tor.rmax + ε || (D >= tor.rmax - ε && dD >= 0)
        incore = tor.rmin > 0 && (D < tor.rmin - ε || (D <= tor.rmin + ε && dD <= 0))
        # the next interval [u, v] of the ray between the surfaces of the tube
        u = outside || incore ? _crossing(tor, point, dir, t, tf, outside ? tor.rmax : tor.rmin, outside) : t
        u > tf && return T(Inf)
        v = min(tf, _crossing(tor, point, dir, u, tf, tor.rmax, false))
        tor.rmin > 0 && (v = min(v, _crossing(tor, point, dir, u, v, tor.rmin, true)))
        # its first part within the ϕ range
        lo, hi = max(u, a1), min(v, b1)
        hi - lo > tol && return lo
        lo, hi = max(u, a2), min(v, b2)
        hi - lo > tol && return lo
        t = v
        t >= tf && return T(Inf)
    end
    return T(Inf)
end

#---Drawing functions-------------------------------------------------------------------------------
#   Grid of nϕ × nθ quads on the outer tube (rmax) and on the inner one (rmin, or collapsed on the
#   central circle if rmin = 0), and the end faces of a ϕ segment. All faces are oriented outwards.
function _torus_grid(tor::Torus, facets)
    nϕ = max(4, round(Int, facets * tor.Δϕ / 2π))
    nθ = max(8, facets ÷ 2)
    return nϕ, nθ
end

function GeometryBasics.coordinates(tor::Torus{T}, facets=36) where T<:AbstractFloat
    nϕ, nθ = _torus_grid(tor, facets)
    points = Point3{T}[]
    for r in (tor.rmax, tor.rmin), i in 0:nϕ, j in 0:nθ
        ϕ, θ = tor.ϕ₀ + tor.Δϕ * i / nϕ, 2π * j / nθ
        ρ = tor.rtor + r * cos(θ)
        push!(points, Point3{T}(ρ * cos(ϕ), ρ * sin(ϕ), r * sin(θ)))
    end
    return points
end

function GeometryBasics.faces(tor::Torus{T}, facets=36) where T<:AbstractFloat
    nϕ, nθ = _torus_grid(tor, facets)
    idx(layer, i, j) = layer * (nϕ + 1) * (nθ + 1) + i * (nθ + 1) + j + 1     # layer 0: rmax, 1: rmin
    quads = QuadFace{Int}[]
    for i in 0:nϕ-1, j in 0:nθ-1
        push!(quads, QuadFace{Int}(idx(0, i, j), idx(0, i+1, j), idx(0, i+1, j+1), idx(0, i, j+1)))
        tor.rmin > 0 && push!(quads, QuadFace{Int}(idx(1, i, j), idx(1, i, j+1), idx(1, i+1, j+1), idx(1, i+1, j)))
    end
    if !tor.fullϕ
        for j in 0:nθ-1
            push!(quads, QuadFace{Int}(idx(0, 0, j), idx(0, 0, j+1), idx(1, 0, j+1), idx(1, 0, j)))      # ϕ₀
            push!(quads, QuadFace{Int}(idx(0, nϕ, j), idx(1, nϕ, j), idx(1, nϕ, j+1), idx(0, nϕ, j+1)))  # ϕ₀ + Δϕ
        end
    end
    return quads
end
