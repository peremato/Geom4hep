#---Trapezoid (TRD)-------------------------------------------------------------
struct Trd{T<:AbstractFloat} <: AbstractShape{T}
    x1::T   # Half-length along x at the surface positioned at -dz
    x2::T   # Half-length along x at the surface positioned at +dz
    y1::T   # Half-length along y at the surface positioned at -dz
    y2::T   # Half-length along y at the surface positioned at +dz
    z::T    # Half-length along z axis
    # cached values
    x2minusx1::T    # Difference between half-legths along x at +dz and -dz
    y2minusy1::T    # Difference between half-legths along y at +dz and -dz
    halfx1plusx2::T # Half-length along x at z = 0
    halfy1plusy2::T # Half-length along y at z = 0
    fx::T           # Tangent of inclination angle along x
    fy::T           # Tangent of inclination angle along y
    calfx::T        # Absolute value of cosine of inclination angle along x
    calfy::T        # Absolute value of cosine of inclination angle along y
    secxz::T        # Reciprocal of fCalfX
    secyz::T        # Reciprocal of fCalfY
    tolerancex::T   # Corrected tolerance for Inside checks on X
    tolerancey::T   # Corrected tolerance for Inside checks on Y

    function Trd{T}(x1, x2, y1, y2, z) where T<:AbstractFloat
        x2minusx1 = x2 - x1
        y2minusy1 = y2 - y1
        halfx1plusx2 = 0.5(x1 + x2)
        halfy1plusy2 = 0.5(y1 + y2)
        fx = 0.5(x1 - x2)/z
        fy = 0.5(y1 - y2)/z
        secxz = sqrt(1. + fx^2)
        secyz = sqrt(1. + fy^2)
        calfx = 1/secxz
        calfy = 1/secyz
        tolerancex = kTolerance(T) * sqrt((x2 - x1)^2 + (2z)^2)
        tolerancey = kTolerance(T) * sqrt((y2 - y1)^2 + (2z)^2)
        new(x1, x2, y1, y2, z, x2minusx1, y2minusy1, halfx1plusx2, halfy1plusy2,
            fx, fy, calfx, calfy, secxz, secyz, tolerancex, tolerancey)
    end
end

#---Basic functions-------------------------------------------------------------
function capacity(trd::Trd{T}) where T<:AbstractFloat
    T(8) * trd.halfx1plusx2 * trd.halfy1plusy2 * trd.z + (T(2)/3) * trd.x2minusx1 * trd.y2minusy1 * trd.z
end
function surface(trd::Trd{T}) where T<:AbstractFloat
    s =  2. * 2. * trd.halfy1plusy2 * sqrt(trd.x2minusx1^2 + (2*trd.z)^2)
    s += 2. * 2. * trd.halfx1plusx2 * sqrt(trd.y2minusy1^2 + (2*trd.z)^2)
    s += 4. * trd.x1 * trd.y1 + 4. * trd.x2 * trd.y2
end
function extent(trd::Trd{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    p = Point3{T}( max(trd.x1, trd.x2), max(trd.y1, trd.y2), trd.z)
    ( -p, p )
end

cross(px, py, vx, vy) = vx * py - vy * px

function inside(trd::Trd{T}, point::Point3{T})::Int64  where T<:AbstractFloat
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

function normal(trd::Trd{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point
    surfaces = 0
    snormal = Vector3{T}(0,0,0)
    distz = abs(z) - trd.z
    xnorm = 1/sqrt(4*trd.z^2 + (trd.x2-trd.x1)^2)
    ynorm = 1/sqrt(4*trd.z^2 + (trd.y2-trd.y1)^2)
    distmx = (-2 * trd.z * x - (trd.x2 - trd.x1) * z - trd.z * (trd.x1 + trd.x2)) * xnorm    
    distpx = ( 2 * trd.z * x - (trd.x2 - trd.x1) * z - trd.z * (trd.x1 + trd.x2)) * xnorm
    distmy = (-2 * trd.z * y - (trd.y2 - trd.y1) * z - trd.z * (trd.y1 + trd.y2)) * ynorm    
    distpy = ( 2 * trd.z * y - (trd.y2 - trd.y1) * z - trd.z * (trd.y1 + trd.y2)) * ynorm
    if abs(distmx) <= kTolerance(T)/2
        surfaces += 1
        snormal += [ -2 * trd.z, 0., -(trd.x2 - trd.x1)] * xnorm
    end
    if abs(distpx) <= kTolerance(T)/2
        surfaces += 1
        snormal += [  2 * trd.z, 0., -(trd.x2 - trd.x1)] * xnorm
    end
    if abs(distmy) <= kTolerance(T)/2
        surfaces += 1
        snormal += [ 0., -2 * trd.z, -(trd.y2 - trd.y1)] * ynorm
    end
    if abs(distpy) <= kTolerance(T)/2
        surfaces += 1
        snormal += [ 0., 2 * trd.z, -(trd.y2 - trd.y1)] * ynorm
    end
    if abs(distz) <= kTolerance(T)/2
        surfaces += 1
        if z > 0.
            snormal += [ 0., 0., 1.]
        else
            snormal += [ 0., 0., -1.]
        end
    end
    if surfaces == 0
        return nothing
    elseif surfaces == 1
        return snormal
    else
        return normalize(snormal)
    end
end

function safety(trd::Trd{T}, point::Point3{T}) where T
    x, y, z = point
    safz = trd.z - abs(z)
    dist = safz
    distx = trd.halfx1plusx2 - trd.fx * z
    safx = (distx - abs(x)) * trd.calfx
    if distx >= 0. && safx < dist
        dist = safx
    end
    disty = trd.halfy1plusy2 - trd.fy * z
    safy = (disty - abs(y)) * trd.calfy
    if disty >= 0. && safy < dist
        dist = safy
    end
    dist
end
safetyToOut(trd::Trd{T}, point::Point3{T}) where T = safety(trd,point)
safetyToIn(trd::Trd{T}, point::Point3{T}) where T = -safety(trd,point)


function faceIntersection(trd::Trd{T}, pos::Point3{T}, dir::Vector3{T}, 
                          forY::Bool, mirroredPoint::Bool, toInside::Bool) where T<:AbstractFloat
    x, y, z = pos
    dx, dy, dz = dir
    if forY
        alongV    = trd.y2minusy1
        v1        = trd.y1
        posV      = y
        posK      = x
        dirV      = dy
        dirK      = dx
        fK        = trd.fx
        fV        = trd.fy
        halfKplus = trd.halfx1plusx2
    else
        alongV    = trd.x2minusx1
        v1        = trd.x1
        posV      = x
        posK      = y
        dirV      = dx
        dirK      = dy
        fK        = trd.fy
        fV        = trd.fx
        halfKplus = trd.halfy1plusy2
    end
    if mirroredPoint
        posV *= -1.
        dirV *= -1.
    end
    ndotv = dirV + fV * dz
    if toInside
        ok = ndotv < 0.
    else
        ok = ndotv > 0.
    end
    if !ok
        return ok, 0.0
    end
    alongZ = 2.0 * trd.z
    # distance from trajectory to face
    dist = (alongZ * (posV - v1) - alongV * (z + trd.z)) / (dz * alongV - dirV * alongZ)
    ok &= dist > kTolerance(T)
    if ok
        # need to make sure z hit falls within bounds
        hitz = z + dist * dz
        ok &= abs(hitz) <= trd.z
        # need to make sure hit on varying dimension falls within bounds
        hitk = posK + dist * dirK
        dK   = halfKplus - fK * hitz; # calculate the width of the varying dimension at hitz
        ok &= abs(hitk) <= dK
        if ok && abs(dist) < kTolerance(T)/2
            dist = 0.
        end
    end
    return ok, dist
end

function distanceToOut(trd::Trd{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat

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

function distanceToIn(trd::Trd{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
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

function Base.show(io::IO, trd::Trd{T}) where T
    print(io, "Trd{$T}",(x1=trd.x1, x2=trd.x2, y1=trd.y1, y2=trd.y2, z=trd.z))
end

function GeometryBasics.coordinates(trd::Trd{T}, facets=6) where {T<:AbstractFloat}
    (; x1, x2, y1, y2, z) = trd
    Point3{T}[(-x1,-y1,-z), ( x1,-y1,-z), (-x1, y1,-z), ( x1, y1,-z),
    (-x2,-y2, z), ( x2,-y2, z), (-x2, y2, z), ( x2, y2, z)]
end

function GeometryBasics.faces(trd::Trd{T}, facets=6) where {T<:AbstractFloat}
    iface = ((1,5,6,2),(3,4,8,7),(1,3,7,5),(2,6,8,4),(1,2,4,3),(5,6,8,7))
    (QuadFace{Int64}(f...) for f in iface)
end


#-------------------------------------------------------------------------------

struct TTrd{T<:AbstractFloat} <: AbstractShape{T}
    x1::T   # Half-length along x at the surface positioned at -dz
    x2::T   # Half-length along x at the surface positioned at +dz
    y1::T   # Half-length along y at the surface positioned at -dz
    y2::T   # Half-length along y at the surface positioned at +dz
    z::T    # Half-length along z axis
    # cached values
    triangles::Vector{Triangle{T}}
    function TTrd{T}(x1, x2, y1, y2, z) where T<:AbstractFloat
        vertices = Point3{T}[(-x1,-y1,-z), ( x1,-y1,-z), (-x1, y1,-z), ( x1, y1,-z),
                             (-x2,-y2, z), ( x2,-y2, z), (-x2, y2, z), ( x2, y2, z)]
        triangles = [Triangle{T}(vertices[i],vertices[j], vertices[k]) for (i,j,k) in box_faces]
        new(x1, x2, y1, y2, z, triangles)
    end
end

function Base.show(io::IO, trd::TTrd{T}) where T
    print(io, "TTrd{$T}",(x1=trd.x1, x2=trd.x2, y1=trd.y1, y2=trd.y2, z=trd.z))
end

function distanceToOut(trd::TTrd{T}, point::Point3{T}, direction::Vector3{T})::T where T<:AbstractFloat
    for triangle in trd.triangles
        dist, ok = intersect(point, direction, triangle)
        ok && return dist
    end
    -1.0
end
#=
function extent(trd::TTrd{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    p = Point3{T}( max(trd.x1, trd.x2), max(trd.y1, trd.y2), trd.z)
    ( -p, p )
end
=#
