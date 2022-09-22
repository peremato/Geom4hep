#Boolean
function distanceToOut_booleanunion(shape, point, dir::Vector3{T})::T where T
    (; left, right, transformation) = shape
    dist = T(0)
    positionA = inside(left, point)
    if positionA != kOutside  # point inside A
        while(true)
            distA = distanceToOut(left, point, dir)
            dist += (distA > 0 && distA < Inf) ? distA : 0
            dist += kPushTolerance(T)
            npoint = point + dist * dir
            lpoint = transformation * npoint
            if inside(right, lpoint) != kOutside # B could be overlapping with A -- and/or connecting A to another part of A
                ldir = transformation * dir
                distB = distanceToOut(right, lpoint, ldir)
                dist += (distB > 0 && distB < Inf) ? distB : 0
                dist += kPushTolerance(T)
                npoint = point + dist * dir
                if inside(left, npoint) == kOutside
                    break
                end
            else
                break
            end
        end
        return dist - kPushTolerance(T)
    else
        lpoint = transformation * point
        positionB = inside(right, lpoint)
        if positionB != kOutside  # point inside B
            ldir = transformation * dir
            while(true)
                distB = distanceToOut(right, lpoint, ldir)
                dist += (distB > 0 && distB < Inf) ? distB : 0
                dist += kPushTolerance(T) # Give a push
                npoint = point + dist * dir
                if inside(left, npoint) != kOutside # A could be overlapping with B -- and/or connecting B to another part of B
                    ldir = transformation * dir
                    distA = distanceToOut(left, npoint, dir)
                    dist += (distA > 0 && distA < Inf) ? distA : 0
                    dist += kPushTolerance(T)
                    npoint = point + dist * dir
                    lpoint = transformation * npoint
                    if inside(right, lpoint) == kOutside
                        break
                    end
                else
                    break
                end
            end
            return dist - kPushTolerance(T)
        else
            return T(-1)
        end
    end 
end
function distanceToOut_booleanintersection(shape, point, dir::Vector3{T})::T where T
    (; left, right, transformation) = shape
    distA = distanceToOut(left, point, dir)
    distB = distanceToOut(right, transformation * point, transformation * dir)
    return min(distA, distB)
end
function distanceToOut_booleansubtraction(shape, point, dir::Vector3{T})::T where T
    (; left, right, transformation) = shape
    distA = distanceToOut(left, point, dir)
    distB = distanceToIn(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

#Plane
function distanceToOut_plane(plane::Plane{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    #The function returns infinity if the plane is not hit from inside, negative
    #if the point is outside
    d = T(Inf)
    ndd = dot(normalize(dir), plane.normal)
    saf = safety(plane, point)
    saf > kTolerance(T) && (d = -T(Inf))
    ndd > 0 && saf < kTolerance(T) && (d = -saf/ndd)
    return d
end

#Trd
function distanceToOut_trd(trd::Trd{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    distance = 0.0

    # hit top Z face?
    trd.calfx 
    trd.calfy 
    safz = trd.z - abs(z)
    out = safz < -kTolerance(T)/2
    distx = trd.halfx1plusx2 - trd.fx * z
    out |= (distx - abs(x)) * trd.calfx < -kTolerance(T)/2
    disty = trd.halfy1plusy2 - trd.fy * z
    out |= (disty - abs(y)) * trd.calfy < -kTolerance(T)/2
    out && return -1.
    if dz > 0.
        distz = (trd.z - z)/abs(dz)
        hitx = abs(x + distz * dx)
        hity = abs(y + distz * dy)
        if hitx <= trd.x2 && hity <= trd.y2
            distance = distz
            abs(distance) < -kTolerance(T)/2 && return 0. 
        end
    end
    # hit bottom Z face?
    if dz < 0.
        distz = (trd.z + z)/abs(dz)
        hitx = abs(x + distz * dx)
        hity = abs(y + distz * dy)
        if hitx <= trd.x1 && hity <= trd.y1
            distance = distz
            abs(distance) < -kTolerance(T)/2 && return 0. 
        end
    end

    # hitting X faces?
    okx, distx = faceIntersection(trd, point, dir, false, false, false)
    if okx
        distance = distx
        return distx < kTolerance(T)/2 ? 0.0 : distx 
    end
    okx, distx = faceIntersection(trd, point, dir, false, true, false)
    if okx
        distance = distx
        return distx < kTolerance(T)/2 ? 0.0 : distx 
    end

    # hitting Y faces?
    oky, disty = faceIntersection(trd, point, dir, true, false, false)
    if oky
        distance = disty
        return disty < kTolerance(T)/2 ? 0.0 : disty 
    end
    oky, disty = faceIntersection(trd, point, dir, true, true, false)
    if oky
        distance = disty
        return disty < kTolerance(T)/2 ? 0.0 : disty 
    end
    return distance 
end
function distanceToOut_trd(trd::TTrd{T}, point::Point3{T}, direction::Vector3{T}) where T<:AbstractFloat
    for triangle in trd.triangles
        dist, ok = intersect(point, direction, triangle)
        ok && return dist
    end
    -1.0
end

#Cone
function distanceToOut_cone(cone::Cone{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
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

#Box
function distanceToOut_box(box::Box{T}, point::Point3{T}, direction::Vector3{T}) where T<:AbstractFloat
    safety = -Inf
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > safety 
            safety = d 
        end
    end
    if safety > kTolerance(T)/2
        return -1.0
    end 
    dist = Inf
    for i in 1:3
        d = (copysign(box.fDimensions[i], direction[i]) - point[i])/direction[i]
        if d < dist
            dist = d
        end
    end
    return dist
    #safety = maximum(abs.(point) - box.fDimensions)
    #distance = minimum((copysign.(box.fDimensions, direction) - point) * (1.0 ./ direction))
    #safety > kTolerance(T)/2 ? -1.0 : distance
end

#Trap
function distanceToOut_trap(trap::Trap{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    (; z, planes) = trap
    dz = dir[3]

    #step 0: if point is outside any plane --> return -1, otherwise initialize at Infinity
    outside = abs(point[3]) > z + kTolerance(T)/2
    outside && return T(-1)
    distance = T(Inf)

    #Step 1: find range of distances along dir between Z-planes (smin, smax)
    dz != 0 && (distance = (copysign(z,dz) - point[3]) / dz)

    #Step 2: find distances for intersections with side planes.
    for i in 1:4
        dist = T(Inf)
        ndir = dot(dir, planes[i].normal)
        safe = safety(planes[i], point)

        safe > kTolerance(T) && return T(-1)
        ndir > 0 && safe < kTolerance(T) && (dist = -safe/ndir)
        dist < distance && (distance = dist)
    end
    return distance
end

#Tube
function distanceToOut_tube(tub::Tube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    w = tub.ϕWedge
    distance::T =  -1.
    done = false

    #First we check all dimensions of the tube, whether points are outside (return -1)
    distz = tub.z - abs(z)
    done |= distz < -kTolerance(T)/2

    rsq  = x*x + y*y
    rdotn = dx*x + dy*y
    crmax = rsq - tub.rmax2

    #if outside of Rmax, return -1
    done |= crmax > kTolerance(T) * tub.rmax
    done && return distance

    if tub.rmin > 0.
        crmin = rsq - tub.rmin2
        done |= crmin < -kTolerance(T) * tub.rmin
        done && return distance
    end
    if tub.Δϕ < 2π
        done |= isOutside(w, x, y)
        done && return distance
    end

    #OK, since we're here, then distance must be non-negative, and the
    #smallest of possible intersections
    !done && (distance = Inf)

    invdirz = 1. / nonzero(dz)
    distz   = dz >= 0 ? (tub.z - z) * invdirz : (-tub.z - z) * invdirz
    !done && dz != 0 && distz < distance && (distance = distz)

    #Find the intersection of the trajectories with the two circles.
    #Here I compute values used in both rmin and rmax calculations.

    invnsq = 1. / nonzero(1. - dz*dz)
    b = invnsq * rdotn

    #rmin
    if tub.rmin > 0.
        crmin *= invnsq
        delta = b * b - crmin
        dist_rmin = -b + (delta > 0. ? -sqrt(delta) : 0.)
        delta > 0. && dist_rmin >= -kTolerance(T)/2 && dist_rmin < distance && (distance = dist_rmin)
    end

    #rmax
    crmax *= invnsq
    delta = b * b - crmax
    dist_rmax = -b + (delta >= 0. ? sqrt(delta) : 0.)
    dist_rmax >= -kTolerance(T)/2 && dist_rmax < distance && (distance = dist_rmax)
 
    #Phi planes

    if tub.Δϕ < 2π
        w = tub.ϕWedge
        p = Point2{T}(x,y)
        d = Vector2{T}(dx,dy)
        if tub.Δϕ < π
            dist_phi = phiPlaneIntersection(p, d, w.along1, w.normal1)
            dist_phi < distance && (distance = dist_phi)
            dist_phi = phiPlaneIntersection(p, d, w.along2, w.normal2)
            dist_phi < distance && (distance = dist_phi)
        elseif tub.Δϕ == π
            dist_phi = phiPlaneIntersection(p, d, w.along2, w.normal2)
            dist_phi < distance && (distance = dist_phi)
        else
            dist_phi = phiPlaneIntersection(p, d, w.along1, w.normal1, ϕgtπ=true)
            dist_phi < distance && (distance = dist_phi)
            dist_phi = phiPlaneIntersection(p, d, w.along2, w.normal2, ϕgtπ=true)
            dist_phi < distance && (distance = dist_phi)
        end
    end
    return distance    
end

#Polycone
function distanceToOut_polycone(pcone::Polycone{T,N}, point::Point3{T}, dir::Vector3{T})::T where {T,N}
    z = point[3]
    dz = dir[3]
    distance = Inf
    if N == 1
        return distanceToOut(pcone.sections[1], point - Vector3{T}(0,0,pcone.zᵢ[1]), dir)
    end
    indexLow  = getSectionIndex(pcone, z - kTolerance(T))
    indexHigh = getSectionIndex(pcone, z + kTolerance(T))
    if indexLow < 0 && indexHigh < 0 
        return -1.
    elseif indexLow < 0 && indexHigh > 0
        index = indexHigh
        inside(pcone.sections[index], point - Vector3{T}(0,0,pcone.zᵢ[index])) == kOutside && return -1.
    elseif indexLow != indexHigh && indexLow > 0
        isin = inside(pcone.sections[indexLow], point - Vector3{T}(0,0,pcone.zᵢ[indexLow]))
        index = isin == kOutside ? indexHigh : indexLow
    else
        index = indexLow
    end

    if index < 0
        return 0.
    else
        inside(pcone.sections[index], point - Vector3{T}(0,0,pcone.zᵢ[index])) == kOutside && return -1.
    end

    distance = T(0)
    increment = dz > 0 ? 1 : -1
    abs(dz) < kTolerance(T) && (increment = 0)
    pn = point
    istep = 0
    while true
        #section surface case
        if distance != 0. || istep < 2
            pn = point + distance * dir - Vector3{T}(0,0,pcone.zᵢ[index])
            inside(pcone.sections[index], pn ) == kOutside && break
        else
            pn -=  Vector3{T}(0,0,pcone.zᵢ[index])
        end
        istep += 1
        dist = distanceToOut(pcone.sections[index], pn, dir)
        dist == -1. && return distance
        if abs(dist) < kTolerance(T)/2
            index1 = index
            if index > 1 && index < N
                index1 += increment
            else
                index == 1 && increment > 0 && (index1 += increment)
                index == N && increment < 0 && (index1 += increment)
            end
            pte = point + (distance + dist) * dir - Vector3{T}(0,0,pcone.zᵢ[index1])
            (inside(pcone.sections[index], pte ) == kOutside || increment == 0) && break
        end
        distance += dist
        index += increment
        (increment == 0 || index < 1 || index > N ) && break
    end
    return distance
end

#CutTube
function distanceToOut_cuttube(ctube::CutTube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    bot, top = ctube.planes
    #Compute distance to cut planes
    distance = min(distanceToOut(bot, point, dir), distanceToOut(top, point, dir))
    #Compute distance to tube
    dtube = distanceToOut(ctube.tube, point, dir)
    dtube < distance && (distance = dtube)
    return distance
end


#Volume
function distanceToOut_volume(agg::Aggregate{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    for pvol in agg.pvolumes
        lpoint = pvol.transformation * point
        inout = inside(pvol.volume.shape, lpoint)
        inout == kSurface && return T(0)
        inout == kInside && return distanceToOut(pvol.volume.shape, lpoint, pvol.transformation * dir)
    end
    return T(-1)
end
