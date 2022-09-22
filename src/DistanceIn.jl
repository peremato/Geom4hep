## Boolean
function distanceToIn_booleanunion(shape::BooleanUnion{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    distA = distanceToIn(left, point, dir)
    distB = distanceToIn(right, transformation * point, transformation * dir)
    return min(distA, distB)
end
# function distanceToIn_booleanintersection(shape::BooleanIntersection{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
#     (; left, right, transformation) = shape

#     positionA = inside(left, point)
#     lpoint = transformation * point
#     positionB = inside(right, lpoint)

#     inLeft = positionA == kInside
#     inRight = positionB == kInside

#     inLeft && inRight && return T(-1)

#     dist = T(0)
#     npoint = point
#     ldir = transformation * dir
#     #  main loop
#     while true
#         d1 = d2 = 0
#         if !inLeft
#             d1 = distanceToIn(left, npoint, dir)
#             d1 = max(d1, kTolerance(T))
#             d1 == T(Inf) && return T(Inf)
#         end
#         if !inRight
#             d2 = distanceToIn(right, lpoint, ldir)
#             d2 = max(d2, kTolerance(T))
#             d2 == T(Inf) && return T(Inf)
#         end
#         if d1 > d2
#             # propagate to left shape
#             dist += d1
#             inleft = true
#             npoint += d1 * dir
#             lpoint = transformation * npoint
#             # check if propagated point is inside right shape, check is done with a little push
#             inRight = inside(right, lpoint + kTolerance(T) * ldir) == kInside
#             inRight && return dist
#             # here inleft=true, inright=false
#         else
#             # propagate to right shape
#             dist += d2
#             inright = true
#             # check if propagated point is inside left shape, check is done with a little push
#             npoint += d2 * dir
#             lpoint = transformation * npoint
#             inLeft = inside(left, npoint + kTolerance(T) * dir) == kInside
#             inLeft && return dist
#         end
#     end
#     return dist
# end
function distanceToIn_booleansubtraction(shape, point::Point3{T}, dir) where T
    (; left, right, transformation) = shape

    lpoint = transformation * point
    ldir = transformation * dir
    positionB = inside(left, lpoint)
    inRight = positionB == kInside

    npoint = point
    dist = T(0)

    while true
        if inRight
            # propagate to outside of '- / RightShape'
            d1 = distanceToOut(right, lpoint, ldir)
            dist += (d1 >= 0 && d1 < Inf) ? d1 + kPushTolerance(T) : 0
            npoint = point + (dist + kPushTolerance(T)) * dir
            lpoint = transformation * npoint
            # now master outside 'B'; check if inside 'A'
            inside(left, npoint) == kInside && distanceToOut(left, npoint, dir) > kPushTolerance(T) && return dist
        end

        # if outside of both we do a max operation master outside '-' and outside '+' ;  find distances to both
        d2 = distanceToIn(left, npoint, dir)
        d2 = max(d2, 0)
        d2 == T(Inf) && return T(Inf)

        d1 = distanceToIn(right, lpoint, ldir)
        if d2 < d1 - kTolerance(T)
          dist += d2 + kPushTolerance(T)
          return dist
        end

        #   propagate to '-'
        dist += (d1 >= 0 && d1 < Inf) ? d1 : 0
        npoint = point + (dist + kPushTolerance(T)) * dir
        lpoint = transformation * npoint
        inRight = true
    end
end

# Plane
function distanceToIn_plane(plane::Plane{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    # The function returns a negative distance for points already inside or
    # direction going outwards (along the normal)
    d = T(-Inf)
    ndd = dot(normalize(dir), plane.normal)
    saf = safety(plane, point)
    ndd < 0 && saf > -kTolerance(T) && (d = -saf/ndd)
    return d
end

# Trd
function distanceToIn_trd(trd::Trd{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    distance::T = Inf

    inz = abs(z) < (trd.z - kTolerance(T)/2)
    distx = trd.halfx1plusx2 - trd.fx * z
    inx = (distx - abs(x)) * trd.calfx > kTolerance(T)/2
    disty = trd.halfy1plusy2 - trd.fy * z
    iny = (disty - abs(y)) * trd.calfy > kTolerance(T)/2
    inside = inx && iny && inz
    if inside
        distance = -1.
    end
    done = inside
    done && return distance
    
    okz = z * dz < 0
    okz &= !inz
    if okz
        distz = (abs(z) - trd.z)/abs(dz)
        hitx = abs(x + distz * dx)
        hity = abs(y + distz * dy)
        okzt = z > (trd.z - kTolerance(T)/2) && hitx <= trd.x2 && hity <= trd.y2
        okzb = z < (-trd.z + kTolerance(T)/2) && hitx <= trd.x1 && hity <= trd.y1
        okz &= ( okzt || okzb )
        if okz
            distance = distz 
        end
    end
    done |= okz
    if done
        return abs(distance) < kTolerance(T)/2 ? 0. : distance
    end

    # hitting X faces?
    okx = false
    if !inx
        okx, distx = faceIntersection(trd, point, dir, false, false, true)
        if okx distance = distx end
        okx, distx = faceIntersection(trd, point, dir, false, true, true)
        if okx distance = distx end
    end
    done |= okx
    if done
        return abs(distance) < kTolerance(T)/2 ? 0. : distance
    end
    if !iny
        oky, disty = faceIntersection(trd, point, dir, true, false, true)
        if oky distance = disty end
        oky, disty = faceIntersection(trd, point, dir, true, true, true)
        if oky distance = disty end
    end
    return abs(distance) < kTolerance(T)/2 ? 0. : distance
end

# Cone
function distanceToIn_cone(cone::Cone{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    (; rmin1, rmax1, rmin2, rmax2, Δϕ, ϕWedge, 
    outerSlope, outerOffset, innerSlope, innerOffset) = cone

    distance = T(Inf)
    done = false

    # outside of Z range and going away?
    distz = abs(z) - cone.z
    done |= (distz > kTolerance(T)/2 && z * dz >= 0) || 
            (abs(distz) < kTolerance(T)/2 && z * dz > 0)
    
    # outside of outer tube and going away?
    rsq  = x * x + y * y
    rmax = rmax1 == rmax2 ? rmax1 : outerOffset + outerSlope * z
    crmax = rsq - rmax * rmax
    rdotn = dx * x + dy * y
    done |= crmax > kTolerance(T) * rmax && rdotn >= 0.
    done && return distance
    
    # Next, check all dimensions of the tube, whether points are inside --> return -1
    distance = T(-1)
    
    # For points inside z-range, return -1
    inside = distz < -kTolerance(T)/2
    inside &= crmax < -kTolerance(T) * rmax
    if rmin1 > 0. || rmin2 > 0.
        rmin = rmin1 == rmin2 ? rmin1 : innerOffset + innerSlope * z
        crmin = rsq - rmin * rmin
        inside &= crmin > kTolerance(T) * rmin
    end
    if Δϕ < 2π
        inside &= isInside(ϕWedge, x, y)
    end
    done |= inside
    done && return distance

    # Next step: check if z-plane is the right entry point (both r,phi
    # should be valid at z-plane crossing)
    distance = T(Inf)
    distz /= nonzero(abs(dz))
    hitx = x + distz * dx
    hity = y + distz * dy
    r2 = (hitx * hitx) + (hity * hity)
    isHittingTopPlane    = z >=  cone.z - kTolerance(T)/2 && r2 <= rmax2 * rmax2 + kTolerance(T)/2
    isHittingBottomPlane = z <= -cone.z + kTolerance(T)/2 && r2 <= rmax1 * rmax1 + kTolerance(T)/2
    okz = isHittingTopPlane || isHittingBottomPlane
    if rmin1 > 0 || rmin2 > 0
        isHittingTopPlane &= r2 >= rmin2 * rmin2 - kTolerance(T)/2
        isHittingBottomPlane &= r2 >= rmin1 * rmin1 - kTolerance(T)/2
        okz &= isHittingTopPlane || isHittingBottomPlane
    end
    if Δϕ < 2π
        okz &= isInside(ϕWedge, hitx, hity)
    end

    !done && okz && (distance = distz)
    done |= okz

    # Next step: intersection of the trajectories with the two conical surfaces
    distOuter = distanceToConicalSurface(cone, point, dir; distToIn=true, innerSurface=false)
    distOuter < distance && (distance = distOuter)
    if rmin1 > 0 || rmin2 > 0
        distInner = distanceToConicalSurface(cone, point, dir; distToIn=true, innerSurface=true)
        distInner < distance && (distance = distInner)
    end
    # Find the intersection of trajectories with Phi planes
    if Δϕ < 2π
        p = Point2{T}(x,y)
        d = Vector2{T}(dx,dy)
        distPhi = phiPlaneIntersection(p, d, ϕWedge.along1, ϕWedge.normal1, toIn=true)
        hit = point + distPhi * dir
        !isnan(distPhi) && abs(z + distPhi * dz) <= cone.z && isInsideR(cone, hit) && isInside(ϕWedge, hit[1],hit[2], tol=kTolerance(T)) && distPhi < distance && (distance = distPhi)  
        if Δϕ != π
            distPhi = phiPlaneIntersection(p, d, ϕWedge.along2, ϕWedge.normal2, toIn=true)
            hit = point + distPhi * dir
            !isnan(distPhi) && abs(z + distPhi * dz) <= cone.z && isInsideR(cone, hit) && isInside(ϕWedge, hit[1],hit[2], tol=kTolerance(T)) && distPhi < distance && (distance = distPhi)
        end
    end
    return distance
end

# Box
function distanceToIn_box(box::Box{T}, point::Point3{T}, direction::Vector3{T}) where T<:AbstractFloat
    distsurf = Inf
    distance = -Inf
    distout = Inf
    for i in 1:3
        din  = (-copysign(box.fDimensions[i],direction[i]) - point[i])/direction[i]
        tout =   copysign(box.fDimensions[i],direction[i]) - point[i]
        dout = tout/direction[i]
        dsur = copysign(tout, direction[i])
        if din > distance 
            distance = din 
        end
        if dout < distout
            distout = dout
        end
        if dsur < distsurf
            distsurf = dsur
        end 
    end
    (distance >= distout || distout <= kTolerance(T)/2 || abs(distsurf) <= kTolerance(T)/2) ? Inf : distance
    #invdir = 1.0 ./ direction
    #tempIn  = -copysign.(box.fDimensions, direction) - point
    #tempOut =  copysign.(box.fDimensions, direction) - point
    #distance = maximum(tempIn * invdir)
    #distout  = minimum(tempOut * invdir)
    #distsurf = abs(minimum(copysign.(tempOut, direction)))
    #(distance >= distout || distout <= kTolerance(T)/2 || distsurf <= kTolerance(T)/2) ? Inf : distance
end

#Tube
function distanceToIn_tube(tub::Tube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    x, y, z = point
    dx, dy, dz = dir
    w = tub.ϕWedge
    distance::T =  Inf
    done = false

    # outside of Z range and going away?
    distz = abs(z) - tub.z
    done |= distz > kTolerance(T)/2 && z * dz >= 0

    # outside of outer tube and going away?
    rsq  = x*x + y*y
    rdotn = dx*x + dy*y
    done |= rsq > tub.tolIrmax2 && rdotn >= 0.
    done && return distance

    # Next, check all dimensions of the tube, whether points are inside --> return -1
    distance = T(-1.)

    # For points inside z-range, return -1
    inside = distz < -kTolerance(T)/2
    inside &= rsq < tub.tolIrmax2
    if tub.rmin > 0.
        inside &= rsq > tub.tolIrmin2
    end
    if tub.Δϕ < 2π
        inside &= isInside(w, x, y, tol=kTolerance(T)/2)
    end
    done |= inside
    done && return distance

    #  Next step: check if z-plane is the right entry point (both r,phi)
    #  should be valid at z-plane crossing)
    distance = Inf
    distz /= nonzero(abs(dz))
    hitx = x + distz * dx
    hity = y + distz * dy
    r2   = hitx * hitx + hity * hity; # radius of intersection with z-plane
    okz  = distz > -kTolerance(T)/2 && z * dz < 0
    okz &= (r2 <= tub.rmax2)
    if tub.rmin > 0.
        okz &= tub.rmin2 <= r2
    end
    if tub.Δϕ < 2π
        okz &= isInside(w, hitx, hity, tol=kTolerance(T)/2)
    end
    !done && okz && (distance = distz)
    done |= okz

    # Next step: is in surface and going inside
    onsurface = rsq  >= tub.tolIrmax2 && rsq <= tub.tolOrmax2 && abs(z) < tub.z + kTolerance(T) && x*dx + y*dy <= 0.
    if tub.rmin > 0.
        onsurface |= rsq  >= tub.tolOrmin2 && rsq <= tub.tolIrmin2 && abs(z) < tub.z + kTolerance(T) && x*dx + y*dy >= 0.
    end
    if tub.Δϕ < 2π
        insector = isInside(w, x, y)
        onsurface &= insector
    end
    !done && onsurface && (distance = 0.)
    done |= onsurface
    done && return distance

    # Next step: intersection of the trajectories with the two circles
    p = Point2{T}(x,y)
    d = Vector2{T}(dx,dy)
    dist_rmax = circleIntersection(p, d, tub.rmax2)
    dist_rmax >= kTolerance(T)/2 && abs(z + dist_rmax * dz) <= tub.z && isInside(w, p + dist_rmax * d) && dist_rmax < distance && (distance = dist_rmax)
    if tub.rmin > 0.
        dist_rmin = circleIntersection(p, d, tub.rmin2, largest=true )
        dist_rmin >= kTolerance(T)/2 && abs(z + dist_rmin * dz) <= tub.z && isInside(w, p + dist_rmin * d) && dist_rmin < distance && (distance = dist_rmin)
    end

    # Calculate intersection between trajectory and the two phi planes
    if tub.Δϕ < 2π
        dist_phi = phiPlaneIntersection(p, d, w.along1, w.normal1, toIn=true)
        hit = p + dist_phi * d
        !isnan(dist_phi) && abs(z + dist_phi * dz) <= tub.z && isInsideR(tub, hit) && isInside(w, hit) && dist_phi < distance && (distance = dist_phi)
        if tub.Δϕ != π
            dist_phi = phiPlaneIntersection(p, d, w.along2, w.normal2, toIn=true)
            hit = p + dist_phi * d
            !isnan(dist_phi) && abs(z + dist_phi * dz) <= tub.z && isInsideR(tub, hit) && isInside(w, hit) && dist_phi < distance && (distance = dist_phi)
        end
    end
    return distance
end

#Volume
function distanceToIn_volume(agg::Aggregate{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    distance = T(Inf)
    for pvol in agg.pvolumes
        lpoint = pvol.transformation * point
        inout = inside(pvol.volume.shape, lpoint)
        inout == kInside && return T(-1)
        ldir = pvol.transformation * dir
        inout == kSurface && normal(pvol.volume.shape, lpoint) ⋅ ldir < 0 && return T(0)
        dist = distanceToIn(pvol.volume.shape, lpoint, ldir)
        dist < distance && (distance = dist)
    end
    return distance
end 

#Polycone
function distanceToIn_polycone(pcone::Polycone{T,N}, point::Point3{T}, dir::Vector3{T})::T where {T,N}
    z = point[3]
    dz = dir[3]
    distance = Inf
    increment = dz > 0 ? 1 : -1
    abs(dz) < kTolerance(T) && (increment = 0)
    index = getSectionIndex(pcone, z)
    index == -1 && (index = 1)
    index == -2 && (index = N)
    while true
        distance = distanceToIn(pcone.sections[index], point - Vector3{T}(0,0,pcone.zᵢ[index]), dir)
        (distance < Inf || increment == 0) && break
        index += increment
        (index < 1 || index > N) && break
    end
    return distance
end

#CutTube
function distanceToIn_cuttube(ctube::CutTube{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
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

#Trap
function distanceToIn_trap(trap::Trap{T}, point::Point3{T}, dir::Vector3{T}) where T<:AbstractFloat
    (; z, planes) = trap
    dz = dir[3]
    z = point[3]

    # step 1.a: input particle is moving away --> return infinity
    signZdir = copysign(1,dz)
    max = signZdir * trap.z - z
    signZdir * max < kTolerance(T)/2 && return Inf

    # step 1.b General case:
    smax = max/dz
    smin = -(signZdir * trap.z + z)/dz 

    # Step 2: find distances for intersections with side planes.
    done = false
    for i in 1:4
        ndir = dot(dir, planes[i].normal)
        safe = safety(planes[i], point)
        dist = -safe/ndir
        
        done |= (safe > kTolerance(T)/2 && ndir >= 0) || (safe > -kTolerance(T) && ndir > 0)

        posPoint = safe > -kTolerance(T)/2
        posDir   = ndir > 0
        !posPoint && posDir && dist < smax && (smax = dist)
        posPoint && !posDir && dist > smin && (smin = dist)
    end
    done && return T(Inf)
    distance = smin <= smax ? smin : Inf
    distance < -kTolerance(T)/2 && return T(-1)
    return distance
end 
