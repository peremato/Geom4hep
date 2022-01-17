#---Cone-------------------------------------------------------------
struct Cone{T<:AbstractFloat} <: AbstractShape{T}
    rmin1::T # inner radius at the surface positioned at -z
    rmax1::T # outer radius at the surface positioned at -z
    rmin2::T # inner radius at the surface positioned at z
    rmax2::T # outer radius at the surface positioned at z
    z::T    # half-length in +z and -z direction
    ϕ₀::T # starting ϕ value (in radians)
    Δϕ::T # delta ϕ value of tube segment (in radians)

    # cached complex values (to avoid recomputation during usage)
    tanRMin::T
    secRMin::T
    tanRMax::T
    secRMax::T 
    innerSlope::T  # "gradient" of inner surface in z direction
    outerSlope::T  # "gradient" of outer surface in z direction
    innerOffset::T
    outerOffset::T
    innerTolerance::T # tolerance on radial direction for inner surface
    outerTolerance::T # tolerance on radial direction for outer surface
    innerConeApex::T
    tanInnerApexAngle::T
    outerConeApex::T
    tanOuterApexAngle::T
    ϕWedge::Wedge{T}
    function Cone{T}(rmin1, rmax1, rmin2, rmax2, z, ϕ₀, Δϕ) where T<:AbstractFloat
        @assert rmax1 >= rmin1  "rmax1 needs to be greater than rmin1"
        @assert rmax2 >= rmin2  "rmax2 needs to be greater than rmin2"
        tanRMin        = (rmin2 - rmin1)/2z
        secRMin        = √(1.0 + tanRMin^2)
        tanRMax        = (rmax2 - rmax1)/2z
        secRMax        = √(1.0 + tanRMax^2)
        innerSlope     = -(rmin1 - rmin2)/2z
        outerSlope     = -(rmax1 - rmax2)/2z
        innerOffset    = rmin2 - innerSlope * z
        outerOffset    = rmax2 - outerSlope * z
        innerTolerance = kTolerance(T) * secRMin
        outerTolerance = kTolerance(T) * secRMax

        if rmin1 == 0. || rmin2 == 0.
            innerConeApex = 0.
            tanInnerApexAngle = rmin1 == 0 ? rmin2 / 2z : rmin1 / 2z
        elseif rmin2 > rmin1
            innerConeApex     = 2z * rmin1 / (rmin2 - rmin1)
            tanInnerApexAngle = rmin2 / (2z + innerConeApex)
        else
            innerConeApex     = 2z * rmin2 / (rmin1 - rmin2)
            tanInnerApexAngle = rmin1 / (2z + innerConeApex)
        end
        if rmax1 == 0. || rmax2 == 0.
            outerConeApex = 0.
            tanOuterApexAngle = rmax1 == 0 ? rmax2 / 2z : rmax1 / 2z
        elseif rmax2 > rmax1
            outerConeApex     = 2z * rmax1 / (rmax2 - rmax1)
            tanOuterApexAngle = rmax2 / (2z + outerConeApex)
        else
            outerConeApex     = 2z * rmax2 / (rmax1 - rmax2)
            tanOuterApexAngle = rmax1 / (2z + outerConeApex)
        end

        new(rmin1, rmax1, rmin2, rmax2, z, ϕ₀, Δϕ, tanRMin, secRMin, tanRMax, secRMax,
            innerSlope, outerSlope, innerOffset, outerOffset, innerTolerance, outerTolerance,
            innerConeApex, tanInnerApexAngle, outerConeApex, tanOuterApexAngle,
            Wedge{T}(ϕ₀,Δϕ))
    end
end

function Base.show(io::IO, cone::Cone{T}) where T
    print(io, "Cone{$T}",(cone.rmin1, cone.rmax1, cone.rmin2, cone.rmax2, cone.z, cone.ϕ₀, cone.Δϕ))
end

#---Basic functions-------------------------------------------------------------
function capacity(cone::Cone{T}) where T<:AbstractFloat
    (; rmin1, rmax1, rmin2, rmax2, z, Δϕ) = cone 
    z * Δϕ / 3. * (rmax1*rmax1 + rmax2*rmax2 + rmax1*rmax2
                 - rmin1*rmin1 - rmin2*rmin2 - rmin1*rmin2)
end

function surface(cone::Cone{T}) where T<:AbstractFloat
    (; rmin1, rmax1, rmin2, rmax2, z, Δϕ) = cone 
    mmin = (rmin1 + rmin2) * 0.5
    mmax = (rmax1 + rmax2) * 0.5
    dmin = (rmin2 - rmin1)
    dmax = (rmax2 - rmax1)
    Δϕ * (mmin * √(dmin * dmin + 4 * z * z) +
          mmax * √(dmax * dmax + 4 * z * z) +
          0.5 * (rmax1 * rmax1 - rmin1 * rmin1 + rmax2 * rmax2 - rmin2 * rmin2))
end

function extent(cone::Cone{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    (; rmin1, rmax1, rmin2, rmax2, z, ϕ₀, Δϕ) = cone 
    extent(Tube{T}(min(rmin1,rmin2), max(rmax1,rmax2),z, ϕ₀, Δϕ))
end

function inside(cone::Cone{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point
    (; rmin1, rmax1, rmin2, rmax2, Δϕ, ϕWedge,
       outerSlope, outerOffset, innerSlope, innerOffset) = cone

    # Check Z
    outside = abs(z) > cone.z + kTolerance(T)/2
    outside && return kOutside
    cinside = abs(z) < cone.z - kTolerance(T)/2

    # Check on RMax
    r2 = x * x + y * y
    rmax = rmax1 == rmax2 ? rmax1 : outerOffset + outerSlope * z
    outside |= r2 > rmax * rmax + kTolerance(T) * rmax
    outside && return kOutside
    cinside &= r2 < rmax * rmax - kTolerance(T) * rmax

    # Check on RMin
    if rmin1 > 0. || rmin2 > 0.
        rmin = rmin1 == rmin2 ? rmin1 : innerOffset + innerSlope * z
        outside |= r2 <= rmin * rmin - kTolerance(T) * rmin
        outside && return kOutside
        cinside &= r2 > rmin * rmin + kTolerance(T) * rmin
    end
    
    # Check on Phi
    if Δϕ < 2π
        outside |= isOutside(ϕWedge, x, y)
        cinside &= isInside(ϕWedge, x, y)
    end
    return outside ? kOutside : cinside ? kInside : kSurface
end

function safetyToIn(cone::Cone{T}, point::Point3{T}) where T<:AbstractFloat
    x,y,z = point
    (; rmin1, rmax1, rmin2, rmax2, Δϕ, ϕWedge, secRMax, secRMin, 
       outerSlope, outerOffset, innerSlope, innerOffset) = cone

    safety = abs(z) - cone.z
    r  = sqrt(x * x + y * y)
    rmax = rmax1 == rmax2 ? rmax1 : outerOffset + outerSlope * z
    safermax = (r - rmax)/secRMax
    safermax > safety && ( safety = safermax )
    
    if rmin1 > 0. || rmin2 > 0.
        rmin = rmin1 == rmin2 ? rmin1 : innerOffset + innerSlope * z
        safermin = (rmin - r)/secRMin
        safermin > safety && ( safety = safermin )
    end
    
    if Δϕ < 2π
        safephi = Δϕ > π ? r : Inf
        dist =  x * ϕWedge.along1[2] - y * ϕWedge.along1[1]
        dist > kTolerance(T)/2 && (safephi = dist)
        dist =  x * ϕWedge.along2[2] - y * ϕWedge.along2[1]
        dist > -kTolerance(T)/2 && dist > safephi && (safephi = dist)
        safephi < Inf && safephi > safety && (safety = safephi)
    end

    return safety
end

function safetyToOut(cone::Cone{T}, point::Point3{T}) where T<:AbstractFloat
    x,y,z = point
    (; rmin1, rmax1, rmin2, rmax2, Δϕ, ϕWedge, secRMax, secRMin, 
    outerSlope, outerOffset, innerSlope, innerOffset) = cone

    safety = cone.z - abs(z)
    r  = sqrt(x * x + y * y)
    rmax = rmax1 == rmax2 ? rmax1 : outerOffset + outerSlope * z
    safermax = (rmax - r)/secRMax
    safermax < safety && ( safety = safermax ) 
    
    if rmin1 > 0. || rmin2 > 0.
        rmin = rmin1 == rmin2 ? rmin1 : innerOffset + innerSlope * z
        safermin = (r - rmin)/secRMin
        safermin < safety && ( safety = safermin )
    end
    
    if Δϕ < 2π
        safephi = safetyToOut(ϕWedge, x, y)
        safephi < safety && (safety = safephi)
    end
    return safety
end

function isOnConicalSurface(cone::Cone{T}, point::Point3{T}; innerSurface::Bool=false)::Bool where T<:AbstractFloat
    x,y,z = point
    r2 = x * x + y * y
    rad = innerSurface ? cone.innerOffset + cone.innerSlope * z : cone.outerOffset + cone.outerSlope * z
    rad2 = rad * rad
    tol = innerSurface ? cone.innerTolerance : cone.outerTolerance
    (r2 >= rad2 - tol * rad) && (r2 <= rad2 + tol * rad) && (abs(z) < cone.z + kTolerance(T)/2)
end

function normal(cone::Cone{T}, point::Point3{T}; innerSurface::Bool=false) where T<:AbstractFloat
    x, y, z = point
    (; rmin1, rmax1, rmin2, rmax2, tanRMin, tanRMax) = cone
    r = √(x^2 + y^2)
    if innerSurface
        rmin1 == rmin2 && rmin1 != 0. ? Vector3{T}(-x, -y, 0.) :  Vector3{T}(-x, -y, r * tanRMin)
    else
        rmax1 == rmax2 && rmax1 != 0. ? Vector3{T}( x,  y, 0.) :  Vector3{T}( x,  y, -r * tanRMax)
    end
end

function distanceToConicalSurface(cone::Cone{T}, point::Point3{T}, dir::Vector3{T}; distToIn::Bool=false, innerSurface::Bool=false)::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    (; rmin1, rmax1, rmin2, rmax2, Δϕ, ϕWedge, secRMax, secRMin, 
    innerConeApex, tanInnerApexAngle, outerConeApex, tanOuterApexAngle) = cone

    distance = Inf
    onConicalSurface = isOnConicalSurface(cone, point, innerSurface=innerSurface)
    if onConicalSurface
        norm = normal(cone, point, innerSurface=innerSurface)
        vdot = dir ⋅ norm
        vdot ≈ 0. && return distance
        if distToIn
            isOnSurfaceAndMovingInside = vdot < 0.
            if Δϕ < 2π
                isOnSurfaceAndMovingInside && isInside(ϕWedge, x, y) && return 0.
            else
                isOnSurfaceAndMovingInside && return 0.
            end
        else
            isOnSurfaceAndMovingOutside = vdot > 0.
            if Δϕ < 2π
                isOnSurfaceAndMovingOutside && isInside(ϕWedge, x, y) && return 0.
            else
                isOnSurfaceAndMovingOutside && return 0.
            end
        end
    end

    ok = false

    pDotV2D = x * dx + y * dy
    if innerSurface
        if rmin1 == rmin2
            a = dx*dx + dy*dy
            b = pDotV2D
            c = x*x + y*y - rmin2*rmin2
        else
            t = tanInnerApexAngle
            newPz = rmin2 > rmin1 ? (z + cone.z + innerConeApex) * t : (z - cone.z - innerConeApex) * t
            a = dx*dx + dy*dy - dz*dz*t*t
            b = pDotV2D - (newPz * dz*t);
            c = x*x + y*y - newPz*newPz
        end
        b2 = b * b
        ac = a * c
        b2 < ac && return NaN
        d2 = b2 - ac
        Δ = √(d2)
        if distToIn
            d2 >= 0. && b >= 0. && (distance = c/(-b - Δ))
            d2 >= 0. && b < 0.  && (distance = (-b + Δ)/a)
        else
            d2 >= 0. && b >= 0. && (distance = (-b - Δ)/a)
            d2 >= 0. && b < 0.  && (distance = c/(-b + Δ))
        end
        distance < 0. && return NaN
        newZ = z + dz * distance
        ok = abs(newZ) < cone.z
    else
        if rmax1 == rmax2
            a = dx*dx + dy*dy
            b = pDotV2D
            c = x*x + y*y - rmax2*rmax2
        else
            t = tanOuterApexAngle
            newPz = rmax2 > rmax1 ? (z + cone.z + outerConeApex) * t : (z - cone.z - outerConeApex) * t
            a = dx*dx + dy*dy - dz*dz*t*t
            b = pDotV2D - (newPz * dz*t);
            c = x*x + y*y - newPz*newPz
        end
        b2 = b * b
        ac = a * c
        b2 < ac && return NaN
        d2 = b2 - ac
        Δ = √(abs(d2))
        if distToIn
            d2 >= 0. && b >= 0. && (distance = (-b - Δ)/a)
            d2 >= 0. && b < 0.  && (distance = c/(-b + Δ))
        else
            d2 >= 0. && b >= 0. && (distance = c/(-b - Δ))
            d2 >= 0. && b < 0.  && (distance = (-b + Δ)/a)
        end
        distance < 0. && return NaN
        newZ = z + dz * distance
        ok = abs(newZ) < cone.z + kTolerance(T)/2
    end
    distance < 0.0 && (distance = Inf)
    if Δϕ < 2π
        hitx = x + distance * dx
        hity = y + distance * dy
        ok &= isInside(ϕWedge, hitx, hity)
    end
    ok ? distance : Inf
end

function distanceToOut(cone::Cone{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    (; rmin1, rmax1, rmin2, rmax2, Δϕ, ϕWedge, secRMax, secRMin, 
    outerSlope, outerOffset, innerSlope, innerOffset) = cone

    distance::T =  -1.
    done = false

    #---First we check all dimensions of the cone, whether points are outside (return -1)
    distz = cone.z - abs(z)
    done |= distz < -kTolerance(T)/2

    r2  = x * x + y * y
    rmax = rmax1 == rmax2 ? rmax1 : outerOffset + outerSlope * z
    crmax = r2 - rmax * rmax

    # if outside of Rmax, return -1
    done |= crmax > kTolerance(T) * rmax
    done && return distance

    if rmin1 > 0. || rmin2 > 0.
        rmin = rmin1 == rmin2 ? rmin1 : innerOffset + innerSlope * z
        crmin = r2 - rmin * rmin
        done |= crmin < -kTolerance(T) * rmin
        done && return distance
    end

    if Δϕ < 2π
        done |= isOutside(ϕWedge, x, y)
        done && return distance
    end

    # OK, since we're here, then distance must be non-negative, and the
    # smallest of possible intersections
    !done && (distance = Inf)

    invdirz = 1. / nonzero(dz)
    distz   = dz >= 0 ? (cone.z - z) * invdirz : (-cone.z - z) * invdirz
    dz != 0 && distz < distance && (distance = distz)

    # Find the intersection of the trajectories with the conical surfaces
    distOuter = distanceToConicalSurface(cone, point, dir; distToIn=false, innerSurface=false)
    distOuter < distance && (distance = distOuter)
    if rmin1 > 0. || rmin2 > 0.
        distInner = distanceToConicalSurface(cone, point, dir; distToIn=false, innerSurface=true)
        distInner < distance && (distance = distInner)
    end

    # Find the intersection of trajectories with Phi planes
    if Δϕ < 2π
        p = Point2{T}(x,y)
        d = Vector2{T}(dx,dy)
        isOnStartPhi = isOnSurface(ϕWedge.along1, ϕWedge.normal1, x, y)
        isOnEndPhi = isOnSurface(ϕWedge.along2, ϕWedge.normal2, x, y)
        cond = (isOnStartPhi && d ⋅ ϕWedge.normal1 < 0.) || (isOnEndPhi && d ⋅ ϕWedge.normal2 < 0.)
        cond && (distance = 0.)
        if Δϕ < π
            dist_phi = phiPlaneIntersection(p, d, ϕWedge.along1, ϕWedge.normal1)
            dist_phi < distance && (distance = dist_phi)
            dist_phi = phiPlaneIntersection(p, d, ϕWedge.along2, ϕWedge.normal2)
            dist_phi < distance && (distance = dist_phi)
        elseif Δϕ == π
            dist_phi = phiPlaneIntersection(p, d, ϕWedge.along2, ϕWedge.normal2)
            dist_phi < distance && (distance = dist_phi)
        else
            dist_phi = phiPlaneIntersection(p, d, ϕWedge.along1, ϕWedge.normal1, ϕgtπ=true)
            dist_phi < distance && (distance = dist_phi)
            dist_phi = phiPlaneIntersection(p, d, ϕWedge.along2, ϕWedge.normal2, ϕgtπ=true)
            dist_phi < distance && (distance = dist_phi)
        end
    end
    return distance    
end


#---Drawing functions------------------------------------------------------------
function GeometryBasics.coordinates(cone::Cone{T}, facets=36) where {T<:AbstractFloat}
    (; rmin1, rmax1, rmin2, rmax2, z, ϕ₀, Δϕ) = cone 
    issector = Δϕ < 2π
    ishollow = rmin1 > 0 || rmin2 > 0
    issector ?  facets =  round(Int64, (facets/2π) * Δϕ) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    nbc = ishollow ? nbv : 1         # Number of centers
    range = 1:(2*nbv + 2*nbc)
    function inner(i)
        return if i <= 2*nbv
            ϕ = T(ϕ₀ + (Δϕ * (((i + 1) ÷ 2) - 1)) / nbf)
            isodd(i) ? Point(rmax2 * cos(ϕ), rmax2 * sin(ϕ), z) : Point(rmax1 * cos(ϕ), rmax1 * sin(ϕ), -z)
        elseif ishollow
            ϕ = T(ϕ₀ + (Δϕ * (((i - 2 * nbv + 1) ÷ 2) - 1)) / nbf)
            isodd(i) ? Point(rmin2 * cos(ϕ), rmin2 * sin(ϕ), z) : Point(rmin1 * cos(ϕ), rmin1 * sin(ϕ), -z)
        elseif i == length(range)
            Point(T(0), T(0), -z)
        elseif i == length(range) - 1
            Point(T(0), T(0), z)
        end
    end
    return (inner(i) for i in range)
  end

  function GeometryBasics.faces(cone::Cone{T}, facets=36) where T<:AbstractFloat
    (; rmin1, rmax1, rmin2, rmax2, z, ϕ₀, Δϕ) = cone 
    issector = Δϕ < 2π
    ishollow = rmin1 > 0 || rmin2 > 0
    issector ?  facets =  round(Int64, (facets/2π) * Δϕ) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    nbc = ishollow ? nbv : 1         # Number of centers

    indexes = Vector{TriangleFace{Int}}()
    for j in 1:nbf
        a,b = 2j-1, 2j
        c,d = !issector && j == nbf ? (1, 2) : (2j+1, 2j+2) 
        push!(indexes, (a,b,d))
        push!(indexes, (d,c,a))
        if ishollow
            a′,b′ = 2j-1+2nbv, 2j+2nbv
            c′,d′ = !issector && j == nbf ? (2nbv+1, 2nbv+2) : (2j+1+2nbv, 2j+2+2nbv)
            # inner wall
            push!(indexes, (a′,d′,b′))
            push!(indexes, (d′,a′,c′))
            # top
            push!(indexes, (a, c ,a′))
            push!(indexes, (c, c′,a′))
            # bottom
            push!(indexes, (b, b′, d))
            push!(indexes, (b′,d′, d))
        else
            a′,b′ = 2nbv+1, 2nbv+2
            # top
            push!(indexes, (a′,a, c))
            # bottom
            push!(indexes, (b′,d, b))
        end
    end
    if issector
        # wedge walls
        a, b, c, d  = ( 1, 2, 2nbv-1, 2nbv)
        a′,b′,c′,d′ = ishollow ? (2nbv+1, 2nbv+2, 4nbv-1, 4nbv ) : (2nbv+1, 2nbv+2, 2nbv+1, 2nbv+2)
        push!(indexes, (a,  b, a′))
        push!(indexes, (b,  b′,a′))
        push!(indexes, (c′, d′,c ))
        push!(indexes, (d′, d, c ))
    end
    return indexes
end
