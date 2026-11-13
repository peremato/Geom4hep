
#---CutTube----------------------------------------------------------------------------------------
struct CutTube{T<:AbstractFloat} <: AbstractShape{T}
    tube::Tube{T}                  # tube (infinite length)
    z::T                           # half z length 
    planes::SVector{2,Plane{T}}    # planes
end

#---Constructor------------------------------------------------------------------------------------
function CutTube{T}(rmin, rmax, z, ϕ₀, Δϕ, bNormal::Vector3{T}, tNormal::Vector3{T}) where T<:AbstractFloat
    CutTube{T}(Tube{T}(rmin, rmax, T(Inf), ϕ₀, Δϕ), z, [Plane{T}(bNormal, Point3{T}(0,0,-z)), Plane{T}(tNormal, Point3{T}(0,0,z))])
end
function CutTube{T}(rmin, rmax, z, ϕ₀, Δϕ, θb, ϕb, θt, ϕt) where T<:AbstractFloat
    CutTube{T}(Tube{T}(rmin, rmax, T(Inf), ϕ₀, Δϕ), z, [Plane{T}(θb, ϕb, Point3{T}(0,0,-z)), Plane{T}( θt, ϕt, Point3{T}(0,0,z))])
end

#---Unilities---------------------------------------------------------------------------------------
z_limit_bottom(ctube::CutTube{T}, r::T, ϕ::T) where T<:AbstractFloat = z_limit(ctube.planes[1], r, ϕ)
z_limit_top(ctube::CutTube{T}, r::T, ϕ::T) where T<:AbstractFloat = z_limit(ctube.planes[2], r, ϕ)

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, ctube::CutTube{T}) where T
    print(io, "CutTube{$T}",(rmin=ctube.tube.rmin, rmax=ctube.tube.rmax, z=ctube.z, ϕ₀=ctube.tube.ϕ₀, Δϕ=ctube.tube.Δϕ, 
               bottom=ctube.planes[1], top=ctube.planes[2]))
end

function GeometryBasics.coordinates(ctube::CutTube{T}, facets=36) where {T<:AbstractFloat}
    tub = ctube.tube
    issector = tub.Δϕ < T(2π)
    ishollow = tub.rmin > 0
    issector ?  facets =  round(Int64, (facets/T(2π)) * tub.Δϕ) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    nbc = ishollow ? nbv : 1         # Number of centers
    z = tub.z
    range = 1:(2*nbv + 2*nbc)
    function inner(i)
        return if i <= 2*nbv
            ϕ = T(tub.ϕ₀ + (tub.Δϕ * (((i + 1) ÷ 2) - 1)) / nbf)
            Point3(tub.rmax * cos(ϕ), tub.rmax * sin(ϕ), z_limit(ctube.planes[i%2+1], tub.rmax, ϕ))
        elseif ishollow
            ϕ = T(tub.ϕ₀ + (tub.Δϕ * (((i - 2 * nbv + 1) ÷ 2) - 1)) / nbf)
            Point3(tub.rmin * cos(ϕ), tub.rmin * sin(ϕ), z_limit(ctube.planes[i%2+1], tub.rmin, ϕ))
        elseif i == length(range)
            Point3(T(0), T(0), -z)
        elseif i == length(range) - 1
            Point3(T(0), T(0), z)
        end
    end
    return (inner(i) for i in range)
end

function GeometryBasics.faces(ctube::CutTube{T}, facets=36) where {T<:AbstractFloat}
    GeometryBasics.faces(ctube.tube, facets)
end

#---Basic functions---------------------------------------------------------------------------------
function extent(ctube::CutTube{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    bot, top = ctube.planes
    z = ctube.z
    (; rmin, rmax, ϕ₀, Δϕ, ϕWedge) = ctube.tube

    dztop = rmax * √(1-top.normal[3]*top.normal[3])/top.normal[3]
    dzbot = -rmax * √(1-bot.normal[3]*bot.normal[3])/bot.normal[3]
    xmin, ymin, zmin = -rmax, -rmax, -z - dzbot
    xmax, ymax, zmax =  rmax,  rmax,  z + dztop
    Δϕ == T(2π) && return (Point3{T}(xmin, ymin, zmin), Point3{T}(xmax, ymax, zmax))

    # The phi cut can reduce the extent in Z
    if !is_inside(ϕWedge, -top.normal[1], -top.normal[2])
        zmax = max(z_limit(top, rmax, ϕ₀), z_limit(top, rmax, ϕ₀ + Δϕ), z_limit(top, rmin,  ϕ₀), z_limit(top, rmin, ϕ₀ + Δϕ))
    end
    if !is_inside(ϕWedge, -bot.normal[1], -bot.normal[2])
        zmin = min(z_limit(bot, rmax, ϕ₀), z_limit(bot, rmax, ϕ₀ + Δϕ), z_limit(bot, rmin,  ϕ₀), z_limit(bot, rmin, ϕ₀ + Δϕ))
    end

    # The phi cut can also reduce the extent in x,y
    c0, c1, s0, s1 = cos(ϕ₀), cos(ϕ₀ + Δϕ), sin(ϕ₀), sin(ϕ₀ + Δϕ)
    is_inside(ϕWedge, T(1), T(0))  || (xmax = max(rmax * c0, rmax * c1, rmin * c0, rmin * c1))
    is_inside(ϕWedge, T(-1), T(0)) || (xmin = min(rmax * c0, rmax * c1, rmin * c0, rmin * c1))
    is_inside(ϕWedge, T(0), T(1))  || (ymax = max(rmax * s0, rmax * s1, rmin * s0, rmin * s1))
    is_inside(ϕWedge, T(0), T(-1)) || (ymin = min(rmax * s0, rmax * s1, rmin * s0, rmin * s1))
    return (Point3{T}(xmin, ymin, zmin), Point3{T}(xmax, ymax, zmax))
end

function inside(ctube::CutTube{T}, point::Point3{T})::Int64  where T<:AbstractFloat
    bot, top = ctube.planes
    tub = ctube.tube
    x, y, z = point

    # Check the cut planes first
    a = safety(bot, point)
    b = safety(top, point)
    outside = a > tolerance(T)/2  || b > tolerance(T)/2 
    outside && return kOutside
    cinside = a < -tolerance(T)/2 && b < -tolerance(T)/2
    inplanes = cinside ? kInside : kSurface
    # Check the tube
    intube = inside(tub, point)
    return inplanes == kSurface && intube != kOutside ? inplanes : intube
end

#   Height between the cut planes at (r, ϕ) is h(r,ϕ) = a + r * (b cosϕ + c sinϕ)
function _height_coefficients(ctube::CutTube{T}) where T<:AbstractFloat
    bot, top = ctube.planes
    a = -top.distance/top.normal[3] + bot.distance/bot.normal[3]
    b = -top.normal[1]/top.normal[3] + bot.normal[1]/bot.normal[3]
    c = -top.normal[2]/top.normal[3] + bot.normal[2]/bot.normal[3]
    return a, b, c
end

function capacity(ctube::CutTube{T}) where T<:AbstractFloat
    (; rmin, rmax, ϕ₀, Δϕ) = ctube.tube
    a, b, c = _height_coefficients(ctube)
    ϕ₁ = ϕ₀ + Δϕ
    a * Δϕ * (rmax^2 - rmin^2) / 2 + (rmax^3 - rmin^3) / 3 * (b * (sin(ϕ₁) - sin(ϕ₀)) + c * (cos(ϕ₀) - cos(ϕ₁)))
end

function surface(ctube::CutTube{T}) where T<:AbstractFloat
    (; rmin, rmax, ϕ₀, Δϕ) = ctube.tube
    bot, top = ctube.planes
    a, b, c = _height_coefficients(ctube)
    ϕ₁ = ϕ₀ + Δϕ
    trig = b * (sin(ϕ₁) - sin(ϕ₀)) + c * (cos(ϕ₀) - cos(ϕ₁))
    lateral(R) = R * (a * Δϕ + R * trig)                      # ∫ h(R,ϕ) R dϕ
    cuts = Δϕ * (rmax^2 - rmin^2) / 2 * (1/abs(bot.normal[3]) + 1/abs(top.normal[3]))
    area = lateral(rmax) + lateral(rmin) + cuts
    if Δϕ < T(2π)
        side(ϕ) = a * (rmax - rmin) + (rmax^2 - rmin^2) / 2 * (b * cos(ϕ) + c * sin(ϕ))   # ∫ h dr
        area += side(ϕ₀) + side(ϕ₁)
    end
    return area
end

function normal(ctube::CutTube{T}, point::Point3{T}) where T<:AbstractFloat
    inside(ctube, point) != kSurface && return zero(Vector3{T})
    n = zero(Vector3{T})
    for plane in ctube.planes
        abs(safety(plane, point)) <= tolerance(T)/2 && (n += plane.normal)
    end
    n += normal(ctube.tube, point)    # lateral surfaces (the tube is infinite in z)
    n ⋅ n > 0 || return zero(Vector3{T})
    return normalize(n)
end

function safety_to_in(ctube::CutTube{T}, point::Point3{T}) where T<:AbstractFloat
    bot, top = ctube.planes
    max(safety_to_in(ctube.tube, point), safety(bot, point), safety(top, point))
end

function safety_to_out(ctube::CutTube{T}, point::Point3{T}) where T<:AbstractFloat
    bot, top = ctube.planes
    min(safety_to_out(ctube.tube, point), -safety(bot, point), -safety(top, point))
end

function distance_to_out(ctube::CutTube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    bot, top = ctube.planes
    # Compute distance to cut planes
    distance = min(distance_to_out(bot, point, dir), distance_to_out(top, point, dir))
    # Compute distance to tube
    dtube = distance_to_out(ctube.tube, point, dir)
    dtube < distance && (distance = dtube)
    return distance
end

function distance_to_in(ctube::CutTube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    bot, top = ctube.planes

    distance = T(Inf)
    # Points already inside have to return negative distance
    pinside = inside(ctube, point)
    pinside == kInside && return T(-1)

    # Outside one of the cut planes and moving away from it: never enters
    for plane in (bot, top)
        safety(plane, point) > tolerance(T)/2 && dot(dir, plane.normal) >= 0 && return T(Inf)
    end

    # Compute distance to cut planes
    d0 = dot(dir, bot.normal) > 0 ? T(-T(Inf)) : distance_to_in(bot, point, dir)
    d1 = dot(dir, top.normal) > 0 ? T(-T(Inf)) : distance_to_in(top, point, dir)
    dplanes = max(d0, d1)
    splanes = max(safety(top,point), safety(bot, point))

    # Mark tracks hitting the planes
    hitplanes = isfinite(dplanes) && splanes > -tolerance(T)
    #!hitplanes && pinside != kInside && return distance

    # Propagate with dplanes only the particles that are hitting
    !hitplanes && (dplanes = 0)
    propagated = point + dplanes * dir

    # Check now which of the propagated points already entering the solid
    intube = inside(ctube.tube, propagated)
    done = hitplanes && intube != kOutside
    done && return abs(dplanes) < tolerance(T) ? zero(T) : dplanes

    # The limit distance for tube crossing cannot exceed the distance to exiting the cut planes
    dexit = min(distance_to_out(bot, propagated, dir), distance_to_out(top, propagated, dir))

    # Compute distance to tube
    dtube = distance_to_in(ctube.tube, propagated, dir)
    dtube < 0 && (dtube = T(0))
    dexit < dtube && (dtube = T(Inf))
    distance = dtube + dplanes
    distance < tolerance(T) && (distance = T(0))
    return distance
end

