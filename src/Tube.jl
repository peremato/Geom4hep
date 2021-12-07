#---Wedge -------------------------------------------------------------
struct Wedge{T<:AbstractFloat} <: AbstractShape{T}
    sphi::T             # starting angle
    dphi::T             # delta angle representing/defining the wedge
    along1::Vector2{T}  # vector along the first plane
    along2::Vector2{T}  # vector along the second plane
    normal1::Vector2{T} # normal vector for first plane
    normal2::Vector2{T} # normal vector for second plane
    function Wedge{T}(sphi, dphi) where T<:AbstractFloat
        along1 = Vector2{T}(cos(sphi), sin(sphi))
        along2 = Vector2{T}(cos(sphi+dphi), sin(sphi+dphi))
        normal1 = Vector2{T}(-sin(sphi), cos(sphi))
        normal2 = Vector2{T}(sin(sphi+dphi), -cos(sphi+dphi))
        new(sphi, dphi, along1, along2, normal1, normal2)
    end
end

function Base.show(io::IO, wed::Wedge{T}) where T
    print(io, "Wedge{$T}",(wed.sphi,wed.dphi))
end

function contains(w::Wedge{T}, p::Point2{T}) where T<:AbstractFloat
    x, y = p
    startx, starty = w.along1
    endx,endy = w.along2
    startCheck = (-x * starty + y * startx)
    endCheck   = (-endx * y + endy * x)
  
    outside = startCheck < 0.
    if (w.dphi < π)
      outside |= endCheck < 0.
    else
      outside &= endCheck < 0.
    end
    return !outside  
end

function isOnSurface(along::Vector2{T}, normal::Vector2{T}, p::Point2{T}) where T<:AbstractFloat  
    along[1]*p[1] + along[2]*p[2] >= 0. && abs(normal[1]*p[1] + normal[2]*p[2]) < kTolerance
end  

function inside(w::Wedge{T}, p::Point2{T}) where T<:AbstractFloat
    startx, starty = w.along1
    endx,endy = w.along2
    startCheck = (-p[1] * w.along1[2] + p[2] * w.along1[1])
    endCheck   = (-w.along2[1] * p[2] + w.along2[2] * p[1])
  
    outside = startCheck < 0.
    if (w.dphi < π)
      outside |= endCheck < 0.
    else
      outside &= endCheck < 0.
    end

    # on right side of half plane && within the right distance to the plane (for both planes) 
    isOnSurface(w.along1, w.normal1, p) || isOnSurface(w.along2, w.normal2, p) ? kSurface : outside ? kOutside : kInside
end

#---Tube -------------------------------------------------------------
struct Tube{T<:AbstractFloat} <: AbstractShape{T}
    rmin::T # inner radius
    rmax::T # outer radius
    z::T    # half-length in +z and -z direction
    sphi::T # starting phi value (in radians)
    dphi::T # delta phi value of tube segment (in radians)
  
    # cached complex values (to avoid recomputation during usage)
    rmin2::T
    rmax2::T
    tolIrmin2::T
    tolOrmin2::T
    tolIrmax2::T
    tolOrmax2::T
    tolIz::T
    tolOz::T
    tolIrmin::T
    tolOrmin::T
    tolIrmax::T
    tolOrmax::T
    maxVal::T
    phiWedge::Wedge{T}
    function Tube{T}(rmin, rmax, z, sphi, dphi) where T<:AbstractFloat
        rmin2 = rmin * rmin
        rmax2 = rmax * rmax
        tolOrmin  = rmin - kTolerance/2
        tolIrmin  = rmin + kTolerance/2
        tolOrmin2 = tolOrmin * tolOrmin
        tolIrmin2 = tolIrmin * tolIrmin
        tolOrmax  = rmax + kTolerance/2
        tolIrmax  = rmax - kTolerance/2
        tolOrmax2 = tolOrmax * tolOrmax
        tolIrmax2 = tolIrmax * tolIrmax
        tolIz  = z - kTolerance/2
        tolOz  = z + kTolerance/2
        maxVal = max(rmax, z)
        new(rmin, rmax, z, sphi, dphi, rmin2, rmax2,
            tolIrmin2, tolOrmin2, tolIrmax2, tolOrmax2, tolIz, tolOz, 
            tolIrmin, tolOrmin, tolIrmax, tolOrmax, maxVal, Wedge{T}(sphi,dphi))
    end
end

function Base.show(io::IO, tub::Tube{T}) where T
    print(io, "Tube{$T}",(tub.rmin, tub.rmax, tub.z, tub.sphi, tub.dphi))
end

#---Basic functions-------------------------------------------------------------
function capacity(tub::Tube{T}) where T<:AbstractFloat
    tub.z * (tub.rmax2 - tub.rmin2) * tub.dphi 
end
function surface(tub::Tube{T}) where T<:AbstractFloat
    toparea = 2 * 0.5 * (tub.rmax2 - tub.rmin2) * tub.dphi
    latphiarea = (tub.dphi < 2π) ? 4. * tub.z * (tub.rmax - tub.rmin) : 0.
    latrinarea = 2. * tub.z * tub.rmin * tub.dphi
    latroutarea = 2. * tub.z * tub.rmax * tub.dphi
    toparea + latphiarea + latrinarea + latroutarea
end
function extent(tub::Tube{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    aMax = [tub.rmax, tub.rmax, tub.z]
    aMin = [-tub.rmax, -tub.rmax, -tub.z]
    if tub.dphi == 2π
        return (Point3{T}(aMin), Point3{T}(aMax))
    end
    # check how many of phi=90, 180, 270, 360deg are outside this tube
    rin  = 0.5 * (tub.rmax + tub.rmin)
    phi0out   = contains(tub.phiWedge, Point2{T}(rin, 0))
    phi90out  = contains(tub.phiWedge, Point2{T}(0, rin))
    phi180out = contains(tub.phiWedge, Point2{T}(-rin, 0))
    phi270out = contains(tub.phiWedge, Point2{T}(0, -rin))
    
    @show phi0out, phi90out, phi180out, phi270out
    # if none of those 4 phis is outside, largest box still required
    if !(phi0out || phi90out || phi180out || phi270out)
        return (Point3{T}(aMin), Point3{T}(aMax))
    end
  
    # some extent(s) of box will be reduced
    # --> think of 4 points A,B,C,D such that A,B are at Rmin, C,D at Rmax
    #     and A,C at startPhi (fSphi), B,D at endPhi (fSphi+fDphi)
    Cx = tub.rmax * cos(tub.sphi)
    Dx = tub.rmax * cos(tub.sphi + tub.dphi)
    Cy = tub.rmax * sin(tub.sphi)
    Dy = tub.rmax * sin(tub.sphi + tub.dphi)
  
    @show Cx, Dx, Cy, Dy
    # then rewrite box sides whenever each one of those phis are not contained in the tube section
    if phi0out; aMax[1] = max(Cx, Dx); end
    if (phi90out); aMax[2] = max(Cy, Dy); end
    if (phi180out) aMin[1] = min(Cx, Dx); end
    if (phi270out) aMin[2] = min(Cy, Dy); end
  
    @show tub.dphi >= π
    if (tub.dphi >= π) 
        return (Point3{T}(aMin), Point3{T}(aMax))
    end
  
    Ax = tub.rmin * cos(tub.sphi)
    Bx = tub.rmin * cos(tub.sphi + tub.dphi)
    Ay = tub.rmin * sin(tub.sphi)
    By = tub.rmin * sin(tub.sphi + tub.dphi)
  
    @show Ax, Bx, Ay, By
    temp     = max(Ax, Bx)
    aMax[1] = temp > aMax[1] ? temp : aMax[1]
  
    temp     = max(Ay, By)
    aMax[2] = temp > aMax[2] ? temp : aMax[2]
  
    temp     = min(Ax, Bx)
    aMin[1] = temp < aMin[1] ? temp : aMin[1]
  
    temp     = min(Ay, By)
    aMin[2] = temp < aMin[2] ? temp : aMin[2]
    return (Point3{T}(aMin), Point3{T}(aMax))
end

function inside(tub::Tube{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point

    # Check Z
    outside = abs(z) > tub.z + kTolerance/2
    outside && return kOutside
    cinside = abs(z) < tub.z - kTolerance/2

    # Check on RMax
    r2 = x * x + y * y
    outside |= r2 > tub.rmax2 + kTolerance * tub.rmax
    outside && return kOutside
    cinside &= r2 < tub.rmax2 - kTolerance * tub.rmax  
    
    # Check on RMin
    if tub.rmin > 0.
        outside |= r2 <= tub.rmin2 - kTolerance * tub.rmin
        outside && return kOutside
        cinside &= r2 > tub.rmin2 + kTolerance * tub.rmin
    end
    
    # Check on Phi
    if tub.dphi < 2π
        insidephi = inside(tub.phiWedge, Point2{T}(x,y))
        outside |= insidephi == kOutside
        cinside  &= insidephi == kInside
    end

    return outside ? kOutside : cinside ? kInside : kSurface
end

function normal(tub::Tube{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point
    norm = zeros(T,3) 
    inside(tub, point) != kSurface && return nothing

    surfaces = 0
    r = sqrt(x*x + y*y)
    inZ = abs(z) - tub.z  <=  kTolerance/2
    inR = r >= (tub.rmin - kTolerance/2) && r <= (tub.rmax + kTolerance/2)
    if inR && abs(z - tub.z) <= kTolerance/2  # top lid, normal along +Z
        norm = [0.,0.,1.]
        surfaces += 1
    end
    if inR && abs(z + tub.z) <= kTolerance/2  # bottom base, normal along -Z
        if surfaces > 0
            norm += [0.,0.,-1]
        else
            norm =  [0.,0.,-1.]
        end
        surfaces += 1
    end

    if tub.rmin > 0.
        if inZ && abs(r - tub.rmin) <= kTolerance/2
            if surfaces > 0
                norm += [-x/r, -y/r, 0.]
            else
                norm =  [-x/r, -y/r, 0.]
            end
            surfaces += 1
        end
    end
    if inZ && abs(r - tub.rmax) <= kTolerance/2
        if surfaces > 0
            norm += [x/r, y/r, 0.]
        else
            norm =  [x/r, y/r, 0.]
        end
        surfaces += 1
    end

    if tub.dphi < 2π
        w = tub.phiWedge
        if inR && isOnSurface(w.along1, w.normal1, Point2{T}(x,y))
            if surfaces > 0
                norm += [-w.normal1[1],-w.normal1[2], 0.]
            else
                norm  = [-w.normal1[1],-w.normal1[2], 0.]
            end
            surfaces += 1
        end

        if inR && isOnSurface(w.along2, w.normal2, Point2{T}(x,y))
            if surfaces > 0
                norm += [-w.normal2[1],-w.normal2[2], 0.]
            else
                norm  = [-w.normal2[1],-w.normal2[2], 0.]
            end
            surfaces += 1
        end
    end
    if surfaces == 0
        return nothing
    elseif surfaces == 1
        return Vector3{T}(norm)
    else
        return Vector3{T}(norm/sqrt(surfaces))
    end
end

function safetyToIn(tub::Tube{T}, point::Point3{T}) where T<:AbstractFloat
    x,y,z = point
    safety = abs(z) - tub.z
    r  = sqrt(x * x + y * y)
    safermax = r - tub.rmax
    safermax > safety && ( safety = safermax ) 
    
    if tub.rmin > 0.
        safermin = tub.rmin - r
        safermin > safety && ( safety = safermin )
    end
    
    if tub.dphi < 2π
        w = tub.phiWedge
        safephi = tub.dphi > π ? r : Inf
        dist =  x * x.along1[2] - y * w.along1[1]
        dist > kTolerance/2 && (safephi = dist)
        dist =  x * w.along2[2] - y * w.along2[1]
        dist > -kTolerance/2 && dist > safephi && (safephi = dist)
        safephi < Inf && safephi > safety && (safety = safephi)
    end

    return safety
end

function safetyToOut(tub::Tube{T}, point::Point3{T}) where T<:AbstractFloat
    x,y,z = point
    safety = tub.z - abs(z)
    r  = sqrt(x * x + y * y)
    safermax = tub.rmax - r
    safermax < safety && ( safety = safermax ) 
    
    if tub.rmin > 0.
        safermin = r - tub.rmin
        safermin < safety && ( safety = safermin )
    end
    
    if tub.dphi < 2π
        w = tub.phiWedge
        safephi = tub.dphi > π ? r : Inf
        dist =  - x * w.along1[2] + y * w.along1[1]
        dist > kTolerance/2 && (safephi = dist)
        dist =  - x * w.along2[2] + y * w.along2[1]
        dist > -kTolerance/2 && dist < safephi && (safephi = dist)
        safephi < safety && (safety = safephi)
    end
    return safety
end


function distanceToOut(tub::Tube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat

    x, y, z = point
    dx, dy, dz = dir
    distance::T =  -1.
    done = false

    #---First we check all dimensions of the tube, whether points are outside (return -1)
    distz = tub.z - abs(z)
    done |= distz < -kTolerance/2

    rsq  = x*x + y*y
    rdotn = dx*x + dy*y
    crmax = rsq - tub.rmax2

    # if outside of Rmax, return -1
    done |= crmax > kTolerance * tub.rmax
    done && return distance

    if tub.rmin > 0.
        crmin = rsq - tub.rmin2
        done |= crmin < -kTolerance * tub.rmin
        done && return distance
    end
    if tub.dphi < 2π
        w = tub.phiWedge
        startCheck = (-x * w.along1[2] + y * w.along1[1])
        endCheck   = (-w.along2[1] * y + w.along2[2] * x)
        done |= startCheck < 0.
        done |= ( tub.dphi < π ? endCheck < 0. : endCheck <= 0.)
        done && return distance
    end

    # OK, since we're here, then distance must be non-negative, and the
    # smallest of possible intersections
    !done && (distance = Inf)

    invdirz = 1. / nonzero(dz)
    distz   = dz >= 0 ? (tub.z - z) * invdirz : (-tub.z - z) * invdirz
    !done && dz != 0 && distz < distance && (distance = distz)

    # Find the intersection of the trajectories with the two circles.
    # Here I compute values used in both rmin and rmax calculations.

    invnsq = 1. / nonzero(1. - dz*dz)
    b = invnsq * rdotn

    # rmin
    if tub.rmin > 0.
        crmin *= invnsq
        delta = b * b - crmin
        dist_rmin = -b + (delta > 0. ? -sqrt(delta) : 0.)
        dist_rmin >= kTolerance/2 && dist_rmin < distance && (distance = dist_rmin)
    end

    # rmax
    crmax *= invnsq
    delta = b * b - crmax
    dist_rmax = -b + (delta >= 0. ? sqrt(delta) : 0.)
    dist_rmax >= kTolerance/2 && dist_rmax < distance && (distance = dist_rmax)
 
    #= Phi planes
    *
    * OK, this is getting weird - the only time I need to
    * check if hit-point falls on the positive direction
    * of the phi-vector is when angle is bigger than PI.
    *
    * Otherwise, any distance I get from there is guaranteed to
    * be larger - so final result would still be correct and no need to
    * check it
    =#

    if tub.dphi < 2π
        w = tub.phiWedge
        if tub.dphi < π
            ok_phi = x * w.normal1[1] + y * w.normal1[2] > 0. 
            dirDotXY = (dy * w.along1[1] - dx * w.along1[2])
            dist_phi = (w.along1[2] * x - w.along1[1] * y ) / nonzero(dirDotXY)
            ok_phi && dist_phi > -kTolerance/2 && dist_phi < distance && (distance = dist_phi)

            ok_phi = x * w.normal2[1] + y * w.normal2[2] > 0. 
            dirDotXY = (dy * w.along2[1] - dx * w.along2[2])
            dist_phi = (w.along2[2] * x - w.along2[1] * y ) / nonzero(dirDotXY)
            ok_phi && dist_phi > -kTolerance/2 && dist_phi < distance && (distance = dist_phi)
        end
    end
    return distance    
end

function distanceToIn(tub::Tube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    distance::T =  Inf
end


function GeometryBasics.coordinates(tub::Tube{T}, facets=36) where {T<:AbstractFloat}
    issector = tub.dphi < 2π
    ishollow = tub.rmin > 0
    issector ?  facets =  round(Int64, (facets/2π) * tub.dphi) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    nbc = ishollow ? nbv : 1         # Number of centers
    z = tub.z
    range = 1:(2*nbv + 2*nbc)
    function inner(i)
        return if i <= 2*nbv
            phi = T((tub.dphi * (((i + 1) ÷ 2) - 1)) / nbf)
            up = ifelse(isodd(i), -z, z)
            Point(tub.rmax * cos(phi), tub.rmax * sin(phi), up)
        elseif ishollow
            phi = T((tub.dphi * (((i - 2 * nbv + 1) ÷ 2) - 1)) / nbf)
            up = ifelse(isodd(i), -z, z)
            Point(tub.rmin * cos(phi), tub.rmin * sin(phi), up)
        elseif i == length(range)
            Point(T(0), T(0), -z)
        elseif i == length(range) - 1
            Point(T(0), T(0), z)
        end
    end
    return (inner(i) for i in range)
  end
  
  function GeometryBasics.faces(tub::Tube{T}, facets=36) where T<:AbstractFloat
    issector = tub.dphi < 2π
    ishollow = tub.rmin > 0
    issector ?  facets =  round(Int64, (facets/2π) * tub.dphi) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    nbc = ishollow ? nbv : 1         # Number of centers
  
    indexes = Vector{TriangleFace{Int}}(undef, nbv * 2 + (ishollow ? nbc * 2 : 0) + (issector ? 2 : 0))
    # External surface
    index = 1
    for j in 1:(nbv - 1)
        indexes[index] = (index + 2, index + 1, index)
        indexes[index + 1] = (index + 3, index + 1, index + 2)
        index += 2
    end
    if issector
        indexes[index]   = (1, 2, index + 2)
        indexes[index+1] = (1, index + 2, index + 3)
        indexes[index+2] = (index + 1, index + 2, index + 3)
        indexes[index+3] = (index + 1, index + 3, index + 2)
        index += 4
    else
        indexes[index] = (1, index + 1, index)
        indexes[index + 1] = (2, index + 1, 1)
        index += 2
    end
    # Internal surface
    if ishollow
        orig = index
        for j in 1:(nbc - 1)
            @show index
            indexes[index] = (index + 2, index + 1, index)
            indexes[index + 1] = (index + 3, index + 1, index + 2)
            index += 2
        end
        if issector
            indexes[index]   = (1, 2, index + 2)
            indexes[index+1] = (1, index + 2, index + 3)
            index += 2
        else
            indexes[index] = (orig, index + 1, index)
            indexes[index + 1] = (orig+1, index + 1, orig)
            index += 2
        end
    end
    # top and bottom
    index = 1
    if ishollow
        for j in 1:(nbv - 1)
            push!(indexes, (index, index + 2 * nbv, index + 2))
            push!(indexes, (index + 2, index + 2 + 2 * nbv, index + 2 * nbv))
            push!(indexes, (index + 1, index + 1 + 2 * nbv, index + 3))
            push!(indexes, (index + 3, index + 3 + 2 * nbv, index + 1 + 2 * nbv))
            index += 2
        end
    else
        for i in 1:(length(indexes) - 2)
            i % 2 == 1 ? push!(indexes, (indexes[i][1], indexes[i][3], 2 * nbv + 2)) :
            push!(indexes, (indexes[i][2], indexes[i][1], 2 * nbv + 1))
        end
    end
    return indexes
  end
  