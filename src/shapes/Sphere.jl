#---Sphere------------------------------------------------------------------------------------------
#   Spherical shell rmin ≤ r ≤ rmax with an optional ϕ segment and θ range, as G4Sphere. The shape is
#   the intersection of up to five regions: the ball r ≤ rmax, the outside of the ball r ≥ rmin (if
#   rmin > 0), the ϕ wedge (if Δϕ < 2π), and the regions θ ≥ θ₀ (if θ₀ > 0) and θ ≤ θ₁ = θ₀ + Δθ (if
#   θ₁ < π). The θ boundaries are cones around the z axis, or the plane z = 0 for θ = π/2. A ray
#   leaves the shape at the first exit from any region, and enters it where it crosses into one
#   region at a point that is inside all the others.
#---------------------------------------------------------------------------------------------------
const THETA_NONE, THETA_PLANE, THETA_CONE = Int32(0), Int32(1), Int32(2)

struct Sphere{T<:AbstractFloat} <: AbstractShape{T}
    rmin::T
    rmax::T
    ϕ₀::T              # starting ϕ (radians)
    Δϕ::T              # ϕ extent (radians)
    θ₀::T              # starting θ (radians)
    Δθ::T              # θ extent (radians)
    cosθ₀::T
    sinθ₀::T
    cosθ₁::T
    sinθ₁::T
    θkind₀::Int32      # boundary at θ₀: THETA_NONE (θ₀ = 0), THETA_PLANE (θ₀ = π/2) or THETA_CONE
    θkind₁::Int32      # boundary at θ₁: THETA_NONE (θ₁ = π), THETA_PLANE (θ₁ = π/2) or THETA_CONE
    fullϕ::Bool
    ϕWedge::Wedge{T}
    function Sphere{T}(rmin, rmax, ϕ₀=0, Δϕ=2π, θ₀=0, Δθ=π) where T<:AbstractFloat
        0 <= rmin < rmax || error("Sphere: invalid radii rmin=$rmin, rmax=$rmax")
        close(a, b) = abs(a - b) < 1e-9
        θ₁ = θ₀ + Δθ
        close(θ₁, π) && (θ₁ = π)
        (0 <= θ₀ < θ₁ <= π) || error("Sphere: invalid θ range ($θ₀, $(θ₀ + Δθ))")
        # cos and sin of the boundary angles, exact for the special values
        cs(θ) = close(θ, 0) ? (1, 0) : close(θ, π/2) ? (0, 1) : close(θ, π) ? (-1, 0) : (cos(θ), sin(θ))
        kind(θ, none) = close(θ, none) ? THETA_NONE : close(θ, π/2) ? THETA_PLANE : THETA_CONE
        (c₀, s₀), (c₁, s₁) = cs(θ₀), cs(θ₁)
        fullϕ = Δϕ >= 2π - 1e-9
        new(rmin, rmax, ϕ₀, fullϕ ? T(2π) : Δϕ, θ₀, θ₁ - θ₀, c₀, s₀, c₁, s₁, kind(θ₀, 0), kind(θ₁, π),
            fullϕ, Wedge{T}(ϕ₀, fullϕ ? T(2π) : Δϕ))
    end
end
Sphere(rmax::T) where T<:AbstractFloat = Sphere{T}(zero(T), rmax)

function Base.show(io::IO, s::Sphere{T}) where T
    print(io, "Sphere{$T}", (s.rmin, s.rmax, s.ϕ₀, s.Δϕ, s.θ₀, s.Δθ))
end

#---Basic functions--------------------------------------------------------------------------------
capacity(s::Sphere) = (s.rmax^3 - s.rmin^3) / 3 * s.Δϕ * (s.cosθ₀ - s.cosθ₁)
function surface(s::Sphere{T}) where T
    area = (s.rmax^2 + s.rmin^2) * s.Δϕ * (s.cosθ₀ - s.cosθ₁)                 # spherical surfaces
    s.fullϕ || (area += (s.rmax^2 - s.rmin^2) * s.Δθ)                         # two ϕ faces
    s.θkind₀ == THETA_NONE || (area += (s.rmax^2 - s.rmin^2) * s.sinθ₀ * s.Δϕ / 2)   # θ cones (or plane)
    s.θkind₁ == THETA_NONE || (area += (s.rmax^2 - s.rmin^2) * s.sinθ₁ * s.Δϕ / 2)
    return area
end
function extent(s::Sphere{T}) where T
    zmax = s.cosθ₀ >= 0 ? s.rmax * s.cosθ₀ : s.rmin * s.cosθ₀
    zmin = s.cosθ₁ <= 0 ? s.rmax * s.cosθ₁ : s.rmin * s.cosθ₁
    r = s.cosθ₀ >= 0 >= s.cosθ₁ ? s.rmax : s.rmax * max(s.sinθ₀, s.sinθ₁)    # θ range with π/2 or not
    (Point3{T}(-r, -r, zmin), Point3{T}(r, r, zmax))
end

#---Signed distances to the boundaries of the regions (positive inside)-----------------------------
#   For the θ regions r⋅sin(θ - θ₀) = ρ cos θ₀ - z sin θ₀ and r⋅sin(θ₁ - θ) = z sin θ₁ - ρ cos θ₁, with
#   ρ = √(x² + y²): exact distances to the cones near them, and with the right sign everywhere.
@inline _rho(p) = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
@inline _rhoxy(p) = sqrt(p[1]^2 + p[2]^2)
@inline _plane1(s::Sphere, p) = s.ϕWedge.normal1[1] * p[1] + s.ϕWedge.normal1[2] * p[2]
@inline _plane2(s::Sphere, p) = s.ϕWedge.normal2[1] * p[1] + s.ϕWedge.normal2[2] * p[2]
@inline function _wedge_distance(s::Sphere{T}, p) where T
    s.fullϕ && return T(Inf)
    d1, d2 = _plane1(s, p), _plane2(s, p)
    s.Δϕ <= π ? min(d1, d2) : max(d1, d2)
end
@inline _theta_distance(c, sn, σ, p) = σ * (_rhoxy(p) * c - p[3] * sn)     # σ = +1: θ ≥ θc, -1: θ ≤ θc
@inline function _sphere_distance(s::Sphere{T}, p) where T
    r = _rho(p)
    d = s.rmax - r
    s.rmin > 0 && (d = min(d, r - s.rmin))
    s.θkind₀ == THETA_NONE || (d = min(d, _theta_distance(s.cosθ₀, s.sinθ₀, 1, p)))
    s.θkind₁ == THETA_NONE || (d = min(d, _theta_distance(s.cosθ₁, s.sinθ₁, -1, p)))
    return min(d, _wedge_distance(s, p))
end

function inside(s::Sphere{T}, point::Point3{T})::Int64 where T<:AbstractFloat
    d = _sphere_distance(s, point)
    d > tolerance(T)/2 ? kInside : d < -tolerance(T)/2 ? kOutside : kSurface
end
safety_to_out(s::Sphere{T}, point::Point3{T}) where T<:AbstractFloat = _sphere_distance(s, point)
safety_to_in(s::Sphere{T}, point::Point3{T}) where T<:AbstractFloat = -_sphere_distance(s, point)

# Outward normal of a θ boundary (the opposite of the gradient of its signed distance)
@inline function _theta_normal(c::T, sn::T, σ, p) where T
    ρ = _rhoxy(p)
    ux, uy = ρ > 0 ? (p[1]/ρ, p[2]/ρ) : (zero(T), zero(T))
    return Vector3{T}(-σ * c * ux, -σ * c * uy, σ * sn)
end

function normal(s::Sphere{T}, point::Point3{T}) where T<:AbstractFloat
    tol = tolerance(T)/2
    r = _rho(point)
    n = zero(Vector3{T})
    r > 0 || return n
    u = Vector3{T}(point[1]/r, point[2]/r, point[3]/r)
    abs(r - s.rmax) < tol && (n += u)
    s.rmin > 0 && abs(r - s.rmin) < tol && (n -= u)
    s.θkind₀ != THETA_NONE && abs(_theta_distance(s.cosθ₀, s.sinθ₀, 1, point)) < tol && (n += _theta_normal(s.cosθ₀, s.sinθ₀, 1, point))
    s.θkind₁ != THETA_NONE && abs(_theta_distance(s.cosθ₁, s.sinθ₁, -1, point)) < tol && (n += _theta_normal(s.cosθ₁, s.sinθ₁, -1, point))
    if !s.fullϕ
        w = s.ϕWedge
        abs(_plane1(s, point)) < tol && w.along1[1]*point[1] + w.along1[2]*point[2] >= 0 && (n -= Vector3{T}(w.normal1[1], w.normal1[2], 0))
        abs(_plane2(s, point)) < tol && w.along2[1]*point[1] + w.along2[2]*point[2] >= 0 && (n -= Vector3{T}(w.normal2[1], w.normal2[2], 0))
    end
    return n == zero(Vector3{T}) ? n : normalize(n)
end

#---Intersections with the surfaces------------------------------------------------------------------
# Roots of |p + t d|² = R² (d is a unit vector), smallest first, NaN if the ray misses the sphere.
# The numerically stable form avoids the cancellation in -b ± √disc.
@inline function _sphere_roots(p, d, R::T) where T
    b = p[1]*d[1] + p[2]*d[2] + p[3]*d[3]
    c = (_rho(p) - R) * (_rho(p) + R)
    disc = b*b - c
    disc < 0 && return T(NaN), T(NaN)
    q = -(b + copysign(sqrt(disc), b))
    q == 0 && return zero(T), zero(T)
    t1, t2 = q, c / q
    return min(t1, t2), max(t1, t2)
end

# Roots of ρ² cos²θc = z² sin²θc (the double cone), NaN if none. Roots on the wrong nappe are
# rejected by the callers. With a direction along a generator (a = 0) the only root is t2.
@inline function _cone_roots(p, d, c::T, sn::T) where T
    c2, s2 = c^2, sn^2
    a = (d[1]^2 + d[2]^2) * c2 - d[3]^2 * s2
    b = (p[1]*d[1] + p[2]*d[2]) * c2 - p[3]*d[3] * s2
    cc = (p[1]^2 + p[2]^2) * c2 - p[3]^2 * s2
    disc = b*b - a*cc
    if disc < 0                                     # rays through the apex give a double root
        disc > -4eps(T) * b*b || return T(NaN), T(NaN)
        disc = zero(T)
    end
    q = -(b + copysign(sqrt(disc), b))
    q == 0 && return T(NaN), T(NaN)
    return q / a, cc / q
end

# Rate of change of the signed distance of a θ region along the direction d, at the point q
@inline function _theta_rate(c::T, sn::T, σ, q, d) where T
    ρ = _rhoxy(q)
    dρ = ρ > 0 ? (q[1]*d[1] + q[2]*d[2]) / ρ : sqrt(d[1]^2 + d[2]^2)
    return σ * (c * dρ - sn * d[3])
end

# A point of the double cone is on the boundary of the θ region if it is on the nappe of θc
@inline _on_nappe(c, q, tol) = c * q[3] >= -tol

# First exit from the θ region (kind, c, sn, σ) at a distance below `dist`
@inline function _theta_exit(dist::T, kind, c::T, sn::T, σ, p, d, tol) where T
    if kind == THETA_PLANE                          # region σ z <= 0, i.e. zside = -σ
        dz = σ * d[3]                               # direction along the outward normal
        dz > 0 && (dist = min(dist, max(zero(T), -σ * p[3]) / dz))
    elseif kind == THETA_CONE
        if abs(_theta_distance(c, sn, σ, p)) < tol  # on the surface: leaving now, or ignore the root at 0
            _theta_rate(c, sn, σ, p, d) < 0 && return zero(T)
        end
        for t in _cone_roots(p, d, c, sn)
            (isnan(t) || t <= tol || t >= dist) && continue
            q = p + t * d
            _on_nappe(c, q, tol) && _theta_rate(c, sn, σ, q, d) < 0 && (dist = t)
        end
    end
    return dist
end

# Closest valid entry into the shape through the boundary of the θ region (kind, c, sn, σ)
@inline function _theta_entry(s::Sphere{T}, best::T, kind, c::T, sn::T, σ, p, d, tol) where T
    if kind == THETA_PLANE
        if σ * d[3] < 0                             # entering the region σ z <= 0
            t = -p[3] / d[3]
            _valid_entry(s, p, d, t, best) && (best = t)
        end
    elseif kind == THETA_CONE
        for t in _cone_roots(p, d, c, sn)
            isnan(t) && continue
            q = p + t * d
            _on_nappe(c, q, tol) && _theta_rate(c, sn, σ, q, d) > 0 && _valid_entry(s, p, d, t, best) && (best = t)
        end
    end
    return best
end

# Crossing of a ϕ plane (normal `nrm` pointing into the wedge, half-plane along `along`): distance t
# and the component of the direction along the inward normal (> 0: entering the wedge)
@inline function _phi_crossing(p, d, along, nrm, ::Type{T}) where T
    dn = nrm[1]*d[1] + nrm[2]*d[2]
    dn == 0 && return T(NaN), dn
    t = -(nrm[1]*p[1] + nrm[2]*p[2]) / dn
    return t, dn
end
@inline _on_halfplane(along, q, tol) = along[1]*q[1] + along[2]*q[2] >= -tol

function distance_to_out(s::Sphere{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    tol = tolerance(T)/2
    _sphere_distance(s, point) < -tol && return T(-1)
    p, d = point, dir
    # outer ball: the largest root
    _, t = _sphere_roots(p, d, s.rmax)
    dist = isnan(t) ? zero(T) : t
    # inner ball: entering it (the smallest root, if ahead)
    if s.rmin > 0
        r = _rho(p)
        pd = p[1]*d[1] + p[2]*d[2] + p[3]*d[3]
        if abs(r - s.rmin) < tol
            pd < 0 && (dist = zero(T))
        else
            t1, _ = _sphere_roots(p, d, s.rmin)
            !isnan(t1) && t1 > 0 && (dist = min(dist, t1))
        end
    end
    # θ boundaries
    dist = _theta_exit(dist, s.θkind₀, s.cosθ₀, s.sinθ₀, 1, p, d, tol)
    dist = _theta_exit(dist, s.θkind₁, s.cosθ₁, s.sinθ₁, -1, p, d, tol)
    # ϕ planes: leaving the wedge through one of its half-planes
    if !s.fullϕ
        w = s.ϕWedge
        for (along, nrm) in ((w.along1, w.normal1), (w.along2, w.normal2))
            t, dn = _phi_crossing(p, d, along, nrm, T)
            (isnan(t) || dn >= 0 || t < -tol) && continue      # not leaving through this plane (ahead)
            t = max(t, zero(T))
            t < dist && _on_halfplane(along, p + t * d, tol) && (dist = t)
        end
    end
    return max(dist, zero(T))
end

# A crossing into one region at distance t is an entry if the point is inside (or on) all the regions
@inline _valid_entry(s::Sphere{T}, p, d, t, best) where T = t >= -tolerance(T) && t < best && _sphere_distance(s, p + t * d) >= -tolerance(T)/2

function distance_to_in(s::Sphere{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    tol = tolerance(T)/2
    _sphere_distance(s, point) > tol && return T(Inf)
    p, d = point, dir
    best = T(Inf)
    # outer sphere, entering the ball: smallest root
    t1, t2 = _sphere_roots(p, d, s.rmax)
    if !isnan(t1)
        _valid_entry(s, p, d, t1, best) && (best = t1)
    else
        return T(Inf)                              # the ray misses the outer sphere
    end
    # inner sphere, leaving the inner ball: largest root
    if s.rmin > 0
        _, u2 = _sphere_roots(p, d, s.rmin)
        !isnan(u2) && _valid_entry(s, p, d, u2, best) && (best = u2)
    end
    # θ boundaries, entering the regions
    best = _theta_entry(s, best, s.θkind₀, s.cosθ₀, s.sinθ₀, 1, p, d, tol)
    best = _theta_entry(s, best, s.θkind₁, s.cosθ₁, s.sinθ₁, -1, p, d, tol)
    # ϕ planes, entering the wedge through one of its half-planes
    if !s.fullϕ
        w = s.ϕWedge
        for (along, nrm) in ((w.along1, w.normal1), (w.along2, w.normal2))
            t, dn = _phi_crossing(p, d, along, nrm, T)
            (isnan(t) || dn <= 0) && continue
            _on_halfplane(along, p + t * d, tol) && _valid_entry(s, p, d, t, best) && (best = t)
        end
    end
    # entering and immediately leaving (grazing an edge): no entry
    best < Inf && distance_to_out(s, p + best * d, d) <= tol && return T(Inf)
    return best
end
