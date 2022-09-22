function inside(shape, point)
    @nospecialize shape
    if shape isa Trap
        inside_trap(shape, point)
    elseif shape isa Trd
        inside_trd(shape, point)
    elseif shape isa Cone
        inside_cone(shape, point)
    elseif shape isa Box
        inside_box(shape, point)
    elseif shape isa Tube
        inside_tube(shape, point)
    elseif shape isa Polycone
        inside_polycone(shape, point)
    elseif shape isa CutTube
        inside_cuttube(shape, point)
    elseif shape isa BooleanUnion
        inside_booleanunion(shape, point)
    elseif shape isa BooleanSubtraction
        inside_booleansubtraction(shape, point)
    elseif shape isa BooleanIntersection
        inside_booleanintersection(shape, point)
    end
end
function inside_booleanunion(shape::BooleanUnion{T, SL, SR}, point::Point3{T})::Int64 where {T,SL,SR}
    (; left, right, transformation) = shape
    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kInside && return kInside

    positionB = inside(right, lpoint)
    positionB == kInside && return kInside

    if positionA == kSurface && positionB == kSurface
        normalA = normal(left, point)
        normalB = normal(right, lpoint) * transformation
        if dot(normalA, normalB) < 0
            return kInside    # touching solids -)(-
        else 
            return kSurface   # overlapping solids =))
        end
    elseif positionA == kSurface || positionB == kSurface
        return kSurface
    else
        return kOutside
    end
end

function inside_booleanintersection(shape::BooleanIntersection{T, SL, SR}, point::Point3{T})::Int64  where {T,SL,SR}
    (; left, right, transformation) = shape
    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kOutside && return kOutside

    positionB = inside(right, lpoint)
    if positionA == kInside && positionB == kInside
        return kInside
    elseif  positionA == kInside && positionB == kSurface ||
            positionA == kSurface && positionB == kInside ||
            positionA == kSurface && positionB == kSurface
        return kSurface
    else
        return kOutside
    end 
end

function inside_booleansubtraction(shape::BooleanSubtraction{T, SL, SR}, point::Point3{T})::Int64  where {T,SL,SR}
    (; left, right, transformation) = shape
    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kOutside && return kOutside

    positionB = inside(right, lpoint)

    if positionA == kInside && positionB == kOutside
        return kInside;
    elseif positionA == kInside && positionB == kSurface ||
           positionB == kOutside && positionA == kSurface
        return kSurface
    else
        return kOutside
    end
end


function inside_trap(trap::Trap{T}, point::Point3{T}) where T<:AbstractFloat
    z = point[3]
    planes = trap.planes
    # inside z?
    outside = abs(z) > trap.z + kTolerance(T)/2
    cinside = abs(z) < trap.z - kTolerance(T)/2

    for i in 1:4
        s = safety(planes[i], point)  #  positive if the point is on the outside halfspace, negative otherwise.
        outside |= s > kTolerance(T)/2
        cinside &= s < -kTolerance(T)/2
    end
    cinside ? kInside : outside ? kOutside : kSurface
end
function inside_cone(cone::Cone{T}, point::Point3{T}) where T<:AbstractFloat
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
function inside_polycone(pcone::Polycone{T,N}, point::Point3{T}) where {T,N}
    x, y, z = point
    indexLow  = getSectionIndex(pcone, z - kTolerance(T))
    indexHigh = getSectionIndex(pcone, z + kTolerance(T))
    indexLow < 0 && indexHigh < 0 && return kOutside
    indexLow < 0 && indexHigh == 1 && return inside(pcone.sections[1], point-Vector3{T}(0,0,pcone.zᵢ[1]))
    indexHigh < 0 && indexLow == N && return inside(pcone.sections[N], point-Vector3{T}(0,0,pcone.zᵢ[N]))
    if indexLow == indexHigh
        return inside(pcone.sections[indexLow], point-Vector3{T}(0,0,pcone.zᵢ[indexLow]))
    else
        insideLow = inside(pcone.sections[indexLow], point-Vector3{T}(0,0,pcone.zᵢ[indexLow]))
        insideHigh = inside(pcone.sections[indexHigh], point-Vector3{T}(0,0,pcone.zᵢ[indexHigh]))
        insideLow == kSurface && insideHigh == kSurface && return kInside
        insideLow == kOutside && insideHigh == kOutside && return kOutside
        return kSurface
    end
end
function inside_box(box::Box{T}, point::Point3{T}) where T<:AbstractFloat
    dist = -Inf
    for i in 1:3
        d = abs(point[i]) - box.fDimensions[i]
        if d  > dist 
            dist = d 
        end
    end
    abs(dist) <= kTolerance(T)/2 ? kSurface : dist < 0.0 ? kInside : kOutside
    #dist = maximum(abs.(point) - box.fDimensions)
    #isapprox(dist, 0.0, atol = kTolerance(T)/2) ? kSurface : dist < 0.0 ? kInside : kOutside
end
function inside_trd(trd::Trd{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point

    # inside z?
    outside = abs(z) > (trd.z + kTolerance(T)/2)
    inside = abs(z) < (trd.z - kTolerance(T)/2)

    # inside x?
    c = cross(abs(x)-trd.x1, z+trd.z, trd.x2-trd.x1, 2*trd.z)
    outside |= (c < -kTolerance(T)/2)
    inside &= (c > kTolerance(T)/2)

    # inside y?
    c = cross(abs(y)-trd.y1, z+trd.z, trd.y2-trd.y1, 2*trd.z)
    outside |= (c < -kTolerance(T)/2)
    inside &= (c > kTolerance(T)/2)

    # return 
    inside ? kInside : outside ? kOutside : kSurface
end
function inside_cuttube(ctube::CutTube{T}, point::Point3{T}) where T<:AbstractFloat
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
function inside_tube(tub::Tube{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point

    # Check Z
    outside = abs(z) > tub.z + kTolerance(T)/2
    outside && return kOutside
    cinside = abs(z) < tub.z - kTolerance(T)/2

    # Check on RMax
    r2 = x * x + y * y
    outside |= r2 > tub.rmax2 + kTolerance(T) * tub.rmax
    outside && return kOutside
    cinside &= r2 < tub.rmax2 - kTolerance(T) * tub.rmax  
    
    # Check on RMin
    if tub.rmin > 0.
        outside |= r2 <= tub.rmin2 - kTolerance(T) * tub.rmin
        outside && return kOutside
        cinside &= r2 > tub.rmin2 + kTolerance(T) * tub.rmin
    end
    
    # Check on Phi
    if tub.Δϕ < 2π
        outside |= isOutside(tub.ϕWedge, x, y)
        cinside &= isInside(tub.ϕWedge, x, y)
    end

    return outside ? kOutside : cinside ? kInside : kSurface
end
