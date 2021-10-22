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
        tolerancex = kTolerance * sqrt((x2 - x1)^2 + (2z)^2)
        tolerancey = kTolerance * sqrt((y2 - y1)^2 + (2z)^2)
        new(x1, x2, y1, y2, z, x2minusx1, y2minusy1, halfx1plusx2, halfy1plusy2,
            fx, fy, calfx, calfy, secxz, secyz, tolerancex, tolerancey)
    end
end

#---Basic functions-------------------------------------------------------------
function capacity(trd::Trd{T}) where T<:AbstractFloat
    8 * trd.halfx1plusx2 * trd.halfy1plusy2 * trd.z + 2//3 * trd.x2minusx1 * trd.y2minusy1 * trd.z
end
function surface(trd::Trd{T}) where T<:AbstractFloat
    s =  2. * 2. * trd.halfy1plusy2 * sqrt(trd.x2minusx1^2 + (2*trd.z)^2)
    s += 2. * 2. * trd.halfx1plusx2 * sqrt(trd.y2minusy1^2 + (2*trd.z)^2)
    s += 4. * trd.x1 * trd.y1 + 4. * trd.x2 * trd.y2
end
function extent(trd::Trd{T}) where T<:AbstractFloat
    p = Point3{T}( max(trd.x1, trd.x2), max(trd.y1, trd.y2), trd.z)
    ( -p, p )
end

cross(px, py, vx, vy) = vx * py - vy * px

function inside(trd::Trd{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point

    # inside z?
    outside = abs(z) > (trd.z + kTolerance/2)
    inside = abs(z) < (trd.z - kTolerance/2)

    # inside x?
    c = cross(abs(x)-trd.x1, z+trd.z, trd.x2-trd.x1, 2*trd.z)
    outside |= (c < -kTolerance/2)
    inside &= (c > kTolerance/2)

    # inside y?
    c = cross(abs(y)-trd.y1, z+trd.z, trd.y2-trd.y1, 2*trd.z)
    outside |= (c < -kTolerance/2)
    inside &= (c > kTolerance/2)
    outside, inside 

    # return 
    inside ? kInside : outside ? kOutside : kSurface
end

function normal(trd::Trd{T}, point::Point3{T}) where T
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
    if abs(distmx) <= kTolerance/2
        surfaces += 1
        snormal += [ -2 * trd.z, 0., -(trd.x2 - trd.x1)] * xnorm
    end
    if abs(distpx) <= kTolerance/2
        surfaces += 1
        snormal += [  2 * trd.z, 0., -(trd.x2 - trd.x1)] * xnorm
    end
    if abs(distmy) <= kTolerance/2
        surfaces += 1
        snormal += [ 0., -2 * trd.z, -(trd.y2 - trd.y1)] * ynorm
    end
    if abs(distpy) <= kTolerance/2
        surfaces += 1
        snormal += [ 0., 2 * trd.z, -(trd.y2 - trd.y1)] * ynorm
    end
    if abs(distz) <= kTolerance/2
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

function distanceToOut(Trd::Trd{T}, point::Point3{T}, direction::Vector3{T}) where T<:AbstractFloat

    x, y, z = point
    dx, dy, dz = direction
    distance::T = 0

    # hit top Z face?
    safz = trd.z - abs(z)
    out = safz < kTolerance/2
    distx = trd.halfx1plusx2 - trd.fx * z
    out |= (distx - abs(x)) * trd.calfx < kTolerance/2
    disty = trd.halfy1plusy2 - trd.fy * z
    out |= (disty - abs(y)) * trd.calfy < kTolerance/2
    out && return -1.
    if dz > 0.
        distz = (trd.z - z)/abs(dz)
        hitx = abs(x + distz * dx)
        hity = abs(y + distz * dy)
        if hitx <= trd.x2 && hity <= trd.y2
            distance = distz
            abs(distance) < kTolerance/2 && return 0. 
        end
    end
    # hit bottom Z face?
    if dz < 0.
        distz = (trd.z + z)/abs(dz)
        hitx = abs(x + distz * dx)
        hity = abs(y + distz * dy)
        if hitx <= trd.x1 && hity <= trd.y1
            distance = distz
            abs(distance) < kTolerance/2 && return 0. 
        end
    end

    # hitting X faces?
 
end
 