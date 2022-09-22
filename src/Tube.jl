#---Tube -------------------------------------------------------------
struct Tube{T<:AbstractFloat} <: AbstractShape{T}
    rmin::T # inner radius
    rmax::T # outer radius
    z::T    # half-length in +z and -z direction
    ϕ₀::T # starting ϕ value (in radians)
    Δϕ::T # delta ϕ value of tube segment (in radians)
  
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
    ϕWedge::Wedge{T}
    function Tube{T}(rmin, rmax, z, ϕ₀, Δϕ) where T<:AbstractFloat
        rmin2 = rmin * rmin
        rmax2 = rmax * rmax
        tolOrmin  = rmin - kTolerance(T)/2
        tolIrmin  = rmin + kTolerance(T)/2
        tolOrmin2 = tolOrmin * tolOrmin
        tolIrmin2 = tolIrmin * tolIrmin
        tolOrmax  = rmax + kTolerance(T)/2
        tolIrmax  = rmax - kTolerance(T)/2
        tolOrmax2 = tolOrmax * tolOrmax
        tolIrmax2 = tolIrmax * tolIrmax
        tolIz  = z - kTolerance(T)/2
        tolOz  = z + kTolerance(T)/2
        maxVal = max(rmax, z)
        new(rmin, rmax, z, ϕ₀, Δϕ, rmin2, rmax2,
            tolIrmin2, tolOrmin2, tolIrmax2, tolOrmax2, tolIz, tolOz, 
            tolIrmin, tolOrmin, tolIrmax, tolOrmax, maxVal, Wedge{T}(ϕ₀,Δϕ))
    end
end
function isInsideR(t::Tube{T}, p::Point2{T}) where T<:AbstractFloat
    r2   = p[1] * p[1] + p[2] * p[2]
    r2 >= t.tolOrmin2 && r2 <= t.tolOrmax2
end

function Base.show(io::IO, tub::Tube{T}) where T
    print(io, "Tube{$T}",(tub.rmin, tub.rmax, tub.z, tub.ϕ₀, tub.Δϕ))
end

#---Basic functions-------------------------------------------------------------
function capacity(tub::Tube{T}) where T<:AbstractFloat
    tub.z * (tub.rmax2 - tub.rmin2) * tub.Δϕ 
end
function surface(tub::Tube{T}) where T<:AbstractFloat
    toparea = 2 * 0.5 * (tub.rmax2 - tub.rmin2) * tub.Δϕ
    latphiarea = (tub.Δϕ < 2π) ? 4. * tub.z * (tub.rmax - tub.rmin) : 0.
    latrinarea = 2. * tub.z * tub.rmin * tub.Δϕ
    latroutarea = 2. * tub.z * tub.rmax * tub.Δϕ
    toparea + latphiarea + latrinarea + latroutarea
end
function extent(tub::Tube{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    aMax = [tub.rmax, tub.rmax, tub.z]
    aMin = [-tub.rmax, -tub.rmax, -tub.z]
    if tub.Δϕ == T(2π)
        return (Point3{T}(aMin), Point3{T}(aMax))
    end
    # check how many of phi=90, 180, 270, 360deg are outside this tube
    rin  = 0.5 * (tub.rmax + tub.rmin)
    ϕ0out   = isOutside(tub.ϕWedge, rin, 0.)
    ϕ90out  = isOutside(tub.ϕWedge, 0., rin)
    ϕ180out = isOutside(tub.ϕWedge, -rin, 0.)
    ϕ270out = isOutside(tub.ϕWedge, 0., -rin)
    
    # if none of those 4 phis is outside, largest box still required
    if !(ϕ0out || ϕ90out || ϕ180out || ϕ270out)
        return (Point3{T}(aMin), Point3{T}(aMax))
    end
  
    # some extent(s) of box will be reduced
    # --> think of 4 points A,B,C,D such that A,B are at Rmin, C,D at Rmax
    #     and A,C at startPhi (fSphi), B,D at endPhi (fSphi+fDphi)
    Cx = tub.rmax * cos(tub.ϕ₀)
    Dx = tub.rmax * cos(tub.ϕ₀ + tub.Δϕ)
    Cy = tub.rmax * sin(tub.ϕ₀)
    Dy = tub.rmax * sin(tub.ϕ₀ + tub.Δϕ)
  
    # then rewrite box sides whenever each one of those phis are not contained in the tube section
    if ϕ0out; aMax[1] = max(Cx, Dx); end
    if (ϕ90out); aMax[2] = max(Cy, Dy); end
    if (ϕ180out) aMin[1] = min(Cx, Dx); end
    if (ϕ270out) aMin[2] = min(Cy, Dy); end
  
    if (tub.Δϕ >= T(π)) 
        return (Point3{T}(aMin), Point3{T}(aMax))
    end
  
    Ax = tub.rmin * cos(tub.ϕ₀)
    Bx = tub.rmin * cos(tub.ϕ₀ + tub.Δϕ)
    Ay = tub.rmin * sin(tub.ϕ₀)
    By = tub.rmin * sin(tub.ϕ₀ + tub.Δϕ)
  
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

function inside(tub::Tube{T}, point::Point3{T})::Int64  where T<:AbstractFloat
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

function normal(tub::Tube{T}, point::Point3{T}) where T<:AbstractFloat
    x, y, z = point
    norm = zeros(T,3) 
    inside(tub, point) != kSurface && return nothing

    surfaces = 0
    r = sqrt(x*x + y*y)
    inZ = abs(z) - tub.z  <=  kTolerance(T)/2
    inR = r >= (tub.rmin - kTolerance(T)/2) && r <= (tub.rmax + kTolerance(T)/2)
    if inR && abs(z - tub.z) <= kTolerance(T)/2  # top lid, normal along +Z
        norm = [0.,0.,1.]
        surfaces += 1
    end
    if inR && abs(z + tub.z) <= kTolerance(T)/2  # bottom base, normal along -Z
        if surfaces > 0
            norm += [0.,0.,-1]
        else
            norm =  [0.,0.,-1.]
        end
        surfaces += 1
    end

    if tub.rmin > 0.
        if inZ && abs(r - tub.rmin) <= kTolerance(T)/2
            if surfaces > 0
                norm += [-x/r, -y/r, 0.]
            else
                norm =  [-x/r, -y/r, 0.]
            end
            surfaces += 1
        end
    end
    if inZ && abs(r - tub.rmax) <= kTolerance(T)/2
        if surfaces > 0
            norm += [x/r, y/r, 0.]
        else
            norm =  [x/r, y/r, 0.]
        end
        surfaces += 1
    end

    if tub.Δϕ < 2π
        w = tub.ϕWedge
        if inR && isOnSurface(w.along1, w.normal1, x, y)
            if surfaces > 0
                norm += [-w.normal1[1],-w.normal1[2], 0.]
            else
                norm  = [-w.normal1[1],-w.normal1[2], 0.]
            end
            surfaces += 1
        end

        if inR && isOnSurface(w.along2, w.normal2, x, y)
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
    
    if tub.Δϕ < 2π
        w = tub.ϕWedge
        safephi = tub.Δϕ > π ? r : Inf
        dist =  x * w.along1[2] - y * w.along1[1]
        dist > kTolerance(T)/2 && (safephi = dist)
        dist =  x * w.along2[2] - y * w.along2[1]
        dist > -kTolerance(T)/2 && dist > safephi && (safephi = dist)
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
    
    if tub.Δϕ < 2π
        w = tub.ϕWedge
        safephi = tub.Δϕ > π ? r : Inf
        dist =  - x * w.along1[2] + y * w.along1[1]
        dist > kTolerance(T)/2 && (safephi = dist)
        dist =  - x * w.along2[2] + y * w.along2[1]
        dist > -kTolerance(T)/2 && dist < safephi && (safephi = dist)
        safephi < safety && (safety = safephi)
    end
    return safety
end

function phiPlaneIntersection(point::Point2{T}, dir::Vector2{T}, along::Vector2{T}, normal::Vector2{T}; ϕgtπ::Bool=false, toIn::Bool=false) where T<:AbstractFloat
    pDotN = point[1] * normal[1] + point[2] * normal[2]
    dDotN = dir[1] * normal[1] + dir[2] * normal[2]
    # point in wrong side
    if (toIn && (pDotN > 0. || pDotN == 0. && dDotN < 0.)) || (!toIn && pDotN <= 0.) 
        return T(NaN)
    end
    dirDotXY = dir[2] * along[1] - dir[1] * along[2]
    dist = (along[2] * point[1] - along[1] * point[2] ) / dirDotXY
    dist < 0. && return T(NaN)
    if ϕgtπ 
        hitx = point[1] + dist * dir[1]
        hity = point[2] + dist * dir[2]
        (hitx * along[1] + hity * along[2]) < 0. && return T(NaN)
    end
    return dist 
end

function circleIntersection(point::Point2{T}, dir::Vector2{T}, R²::T; largest::Bool=false) where T<:AbstractFloat
    # returns the distance to the circle (in units of norm(dir))
    r2 = point[1] * point[1] + point[2] * point[2]
    dr2 = dir[1] * dir[1] + dir[2] * dir[2]
    b = (dir[1] * point[1] + dir[2] * point[2])/dr2
    c = (r2 - R²)/dr2
    Δ = b * b - c
    Δ < 0. && return T(NaN)
    dist = largest ? -b + sqrt(Δ) : -b - sqrt(Δ)
    dist < 0. && return T(NaN)
    return dist
end


function GeometryBasics.coordinates(tub::Tube{T}, facets=36) where {T<:AbstractFloat}
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
            up = ifelse(isodd(i), z, -z)
            Point(tub.rmax * cos(ϕ), tub.rmax * sin(ϕ), up)
        elseif ishollow
            ϕ = T(tub.ϕ₀ + (tub.Δϕ * (((i - 2 * nbv + 1) ÷ 2) - 1)) / nbf)
            up = ifelse(isodd(i), z, -z)
            Point(tub.rmin * cos(ϕ), tub.rmin * sin(ϕ), up)
        elseif i == length(range)
            Point(T(0), T(0), -z)
        elseif i == length(range) - 1
            Point(T(0), T(0), z)
        end
    end
    return (inner(i) for i in range)
  end
  
  function GeometryBasics.faces(tub::Tube{T}, facets=36) where T<:AbstractFloat
    issector = tub.Δϕ < 2π
    ishollow = tub.rmin > 0
    issector ?  facets =  round(Int64, (facets/2π) * tub.Δϕ) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    nbc = ishollow ? nbv : 1         # Number of centers
    indexes = Vector{QuadFace{Int64}}()
    for j in 1:nbf
        a,b = 2j-1, 2j
        c,d = !issector && j == nbf ? (1, 2) : (2j+1, 2j+2) 
        push!(indexes, (a,b,d,c))
        if ishollow
            a′,b′ = 2j-1+2nbv, 2j+2nbv
            c′,d′ = !issector && j == nbf ? (2nbv+1, 2nbv+2) : (2j+1+2nbv, 2j+2+2nbv)
            # inner wall
            push!(indexes, (a′,b′,d′,c′))
            # top
            push!(indexes, (c, c′, a′, a))
            # bottom
            push!(indexes, (b, b′, d′, d))
        else
            a′,b′ = 2nbv+1, 2nbv+2
            # top
            push!(indexes, (a′,a, c, c))
            # bottom
            push!(indexes, (b′,d, b, b))
        end
    end
    if issector
        # wedge walls
        a, b, c, d  = ( 1, 2, 2nbv-1, 2nbv)
        a′,b′,c′,d′ = ishollow ? (2nbv+1, 2nbv+2, 4nbv-1, 4nbv ) : (2nbv+1, 2nbv+2, 2nbv+1, 2nbv+2)
        push!(indexes, (a,  b, b′, a′))
        push!(indexes, (c′, d′, d, c ))
    end
    return indexes
end

function GeometryBasics.normals(tub::Tube{T}, facets=36) where {T<:AbstractFloat}
    issector = tub.Δϕ < 2π
    ishollow = tub.rmin > 0
    issector ?  facets =  round(Int64, (facets/2π) * tub.Δϕ) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    nbc = ishollow ? nbv : 1         # Number of centers
    range = 1:(2*nbv + 2*nbc)
    function inner(i)
        return if i <= 2*nbv
            ϕ = T(tub.ϕ₀ + (tub.Δϕ * (((i + 1) ÷ 2) - 1)) / nbf)
            up = ifelse(isodd(i), 1/√2, -1/√2)
            Vector3(cos(ϕ)/√2, sin(ϕ)/√2, up)
        elseif ishollow
            ϕ = T(tub.ϕ₀ + (tub.Δϕ * (((i + 1) ÷ 2) - 1)) / nbf)
            up = ifelse(isodd(i), 1/√2, -1/√2)
            Vector3(-cos(ϕ)/√2, -sin(ϕ)/√2, up)
        elseif i == length(range)
            Vector3(T(0),T(0), -1.)
        elseif i == length(range) - 1
            Vector3(T(0), T(0), 1.)
        end
    end
    return (inner(i) for i in range)
end
