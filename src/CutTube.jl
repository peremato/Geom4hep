include("Plane.jl")

#---CutTube----------------------------------------------------------------------------------------
struct CutTube{T<:AbstractFloat} <: AbstractShape{T}
    tube::Tube{T}                  # tube (infinite length)
    z::T                           # half z length 
    planes::SVector{2,Plane{T}}    # planes
end

#---Constructor------------------------------------------------------------------------------------
function CutTube{T}(rmin, rmax, z, ϕ₀, Δϕ, bNormal::Vector3{T}, tNormal::Vector3{T}) where T<:AbstractFloat
    CutTube{T}(Tube{T}(rmin, rmax, Inf, ϕ₀, Δϕ), z, [Plane{T}(bNormal, Point3{T}(0,0,-z)), Plane{T}(tNormal, Point3{T}(0,0,z))])
end
function CutTube{T}(rmin, rmax, z, ϕ₀, Δϕ, θb, ϕb, θt, ϕt) where T<:AbstractFloat
    CutTube{T}(Tube{T}(rmin, rmax, Inf, ϕ₀, Δϕ), z, [Plane{T}(θb, ϕb, Point3{T}(0,0,-z)), Plane{T}( θt, ϕt, Point3{T}(0,0,z))])
end

#---Unilities---------------------------------------------------------------------------------------
zLimitBottom(ctube::CutTube{T}, r::T, ϕ::T) where T<:AbstractFloat = zLimit(ctube.planes[1], r, ϕ)
zLimitTop(ctube::CutTube{T}, r::T, ϕ::T) where T<:AbstractFloat = zLimit(ctube.planes[2], r, ϕ)

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, ctube::CutTube{T}) where T
    print(io, "CutTube{$T}",(rmin=ctube.tube.rmin, rmax=ctube.tube.rmax, z=ctube.z, ϕ₀=ctube.tube.ϕ₀, Δϕ=ctube.tube.Δϕ, 
               bottom=ctube.planes[1], top=ctube.planes[2]))
end

function GeometryBasics.coordinates(ctube::CutTube{T}, facets=36) where {T<:AbstractFloat}
    tub = ctube.tube
    issector = tub.Δϕ < 2π
    ishollow = tub.rmin > 0
    issector ?  facets =  round(Int64, (facets/2π) * tub.Δϕ) : nothing
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
            Point(tub.rmax * cos(ϕ), tub.rmax * sin(ϕ), zLimit(ctube.planes[i%2+1], tub.rmax, ϕ))
        elseif ishollow
            ϕ = T(tub.ϕ₀ + (tub.Δϕ * (((i - 2 * nbv + 1) ÷ 2) - 1)) / nbf)
            Point(tub.rmin * cos(ϕ), tub.rmin * sin(ϕ), zLimit(ctube.planes[i%2+1], tub.rmin, ϕ))
        elseif i == length(range)
            Point(T(0), T(0), -z)
        elseif i == length(range) - 1
            Point(T(0), T(0), z)
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
    aMax = [ rmax,  rmax, z + dztop]
    aMin = [-rmax, -rmax, -z - dzbot]
    if Δϕ == T(2π)
        return (Point3{T}(aMin), Point3{T}(aMax))
    end
    
    # The phi cut can reduce the extent in Z
    if !isInside(ϕWedge, -top.normal[1], -top.normal[2])
        aMax[3] = max(zLimit(top, rmax, ϕ₀), zLimit(top, rmax, ϕ₀ + Δϕ), zLimit(top, rmin,  ϕ₀), zLimit(top, rmin, ϕ₀ + Δϕ)) 
    end
    if !isInside(ϕWedge, -bot.normal[1], -bot.normal[2])
        aMin[3] = min(zLimit(bot, rmax, ϕ₀), zLimit(bot, rmax, ϕ₀ + Δϕ), zLimit(bot, rmin,  ϕ₀), zLimit(bot, rmin, ϕ₀ + Δϕ)) 
    end

    # The phi cut can also reduce the extent in x,y
    if !isInside(ϕWedge, T(1), T(0))
        aMax[1] = max(rmax * cos(ϕ₀), rmax * cos(ϕ₀ + Δϕ), rmin * cos(ϕ₀), rmin * cos(ϕ₀ + Δϕ))
    end
    if !isInside(ϕWedge, T(-1), T(0))
        aMin[1] = min(rmax * cos(ϕ₀), rmax * cos(ϕ₀ + Δϕ), rmin * cos(ϕ₀), rmin * cos(ϕ₀ + Δϕ))
    end
    if !isInside(ϕWedge, T(0), T(1))
        aMax[2] = max(rmax * cos(ϕ₀), rmax * cos(ϕ₀ + Δϕ), rmin * cos(ϕ₀), rmin * cos(ϕ₀ + Δϕ))
    end
    if !isInside(ϕWedge, T(0), T(-1))
        aMin[2] = min(rmax * cos(ϕ₀), rmax * cos(ϕ₀ + Δϕ), rmin * cos(ϕ₀), rmin * cos(ϕ₀ + Δϕ))
    end
    return (Point3{T}(aMin), Point3{T}(aMax))
end

@override function inside(ctube::CutTube{T}, point::Point3{T}) where T<:AbstractFloat
    bot, top = ctube.planes
    tub = ctube.tube
    x, y, z = point

    # Check the cut planes first
    a = safety(bot, point)
    b = safety(top, point)
    outside = a > kTolerance(T)/2  || b > kTolerance(T)/2 
    outside && return kOutside
    cinside = a < -kTolerance(T)/2 && b < -kTolerance(T)/2
    inplanes = cinside ? kInside : kSurface
    # Check the tube
    intube = inside(tub, point)
    return inplanes == kSurface && intube != kOutside ? inplanes : intube
end

@override function safetyToIn(ctube::CutTube{T}, point::Point3{T}) where T<:AbstractFloat

end

@override function distanceToOut(ctube::CutTube{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    bot, top = ctube.planes
    # Compute distance to cut planes
    distance = min(distanceToOut(bot, point, dir), distanceToOut(top, point, dir))
    # Compute distance to tube
    dtube = distanceToOut(ctube.tube, point, dir)
    dtube < distance && (distance = dtube)
    return distance
end

@override function distanceToIn(ctube::CutTube{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    bot, top = ctube.planes

    distance = T(Inf)
    # Points already inside have to return negative distance
    pinside = inside(ctube, point)
    pinside == kInside && return T(-1)

    # Compute distance to cut planes
    d0 = dot(dir, bot.normal) > 0 ? T(-Inf) : distanceToIn(bot, point, dir)
    d1 = dot(dir, top.normal) > 0 ? T(-Inf) : distanceToIn(top, point, dir)
    dplanes = max(d0, d1)
    splanes = max(safety(top,point), safety(bot, point))

    # Mark tracks hitting the planes
    hitplanes = dplanes < T(Inf) && splanes > -kTolerance(T)
    #!hitplanes && pinside != kInside && return distance

    # Propagate with dplanes only the particles that are hitting
    !hitplanes && (dplanes = 0)
    propagated = point + dplanes * dir

    # Check now which of the propagated points already entering the solid
    intube = inside(ctube.tube, propagated)
    done = hitplanes && intube != kOutside
    done && return abs(dplanes) < kTolerance(T) ? 0 : dplanes

    # The limit distance for tube crossing cannot exceed the distance to exiting the cut planes
    dexit = min(distanceToOut(bot, propagated, dir), distanceToOut(top, propagated, dir))

    # Compute distance to tube
    dtube = distanceToIn(ctube.tube, propagated, dir)
    dtube < 0 && (dtube = T(0))
    dexit < dtube && (dtube = T(Inf))
    distance = dtube + dplanes
    distance < kTolerance(T) && (distance = T(0))
    return distance
end

