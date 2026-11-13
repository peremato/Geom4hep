#---Polycone----------------------------------------------------------------------------------------
#   A sequence of cone sections along z. The storage of the sections is a type parameter: an SVector
#   for shapes built by the user, or a view into the shared arrays of a flat GeometryModel.
struct Polycone{T<:AbstractFloat, VC<:AbstractVector{Cone{T}}, VZ<:AbstractVector{T}} <: AbstractShape{T}
    sections::VC                   # cone sections (non-zero length)
    zᵢ::VZ                         # z offset of the centre of each section
    ϕ₀::T                          # starting ϕ value (in radians)
    Δϕ::T                          # delta ϕ value of tube segment (in radians)
end

#---Constructor------------------------------------------------------------------------------------
function Polycone{T}(rmin::Vector, rmax::Vector, z::Vector, ϕ₀, Δϕ) where T<:AbstractFloat
    @assert length(rmin) == length(rmax) == length(z) "Vectors for rmin, rmax and z are different size"
    L = length(rmin)
    sections = [Cone{T}(rmin[i],rmax[i], rmin[i+1], rmax[i+1], (z[i+1]-z[i])/2, ϕ₀, Δϕ) for i in 1:L-1 if (z[i+1] - z[i]) > tolerance(T)]
    zoffs = [T((z[i]+z[i+1])/2) for i in 1:L-1 if (z[i+1] - z[i]) > tolerance(T)]
    N = length(sections)
    Polycone{T}(SVector{N,Cone{T}}(sections), SVector{N,T}(zoffs), T(ϕ₀), T(Δϕ))
end
Polycone{T}(sections::VC, zᵢ::VZ, ϕ₀, Δϕ) where {T, VC<:AbstractVector{Cone{T}}, VZ<:AbstractVector{T}} =
    Polycone{T,VC,VZ}(sections, zᵢ, ϕ₀, Δϕ)

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, pcone::Polycone{T}) where T
    (; sections, zᵢ) = pcone
    z = [[zᵢ[i] - sections[i].z for i in eachindex(sections)]; zᵢ[end] + sections[end].z]
    print(io, "Polycone{$T}",(nsections=length(sections), z=z, ϕ₀=pcone.ϕ₀, Δϕ=pcone.Δϕ))
end


function GeometryBasics.mesh(pcone::Polycone{T}; facets::Integer=64) where T
    N = length(pcone.sections)
    merge([Mesh(map(x-> _gbpoint(x + Vector3{T}(0,0,pcone.zᵢ[i])), coordinates(pcone.sections[i], facets)), collect(faces(pcone.sections[i], facets))) for i in 1:N])
end

#---Basic functionality-----------------------------------------------------------------------------
nz(pcone::Polycone) = length(pcone.sections) + 1
nsections(pcone::Polycone) = length(pcone.sections)
function section_index(pcone::Polycone{T}, zpos::T) where T
    N = length(pcone.sections)
    (; sections, zᵢ) = pcone
    zpos < zᵢ[1] - sections[1].z && return -1
    for i  in 1:N
        Δz = sections[i].z
        z  = zᵢ[i]
        zpos >= z - Δz && zpos <= z + Δz && return i
    end
    return -2
end

function capacity(pcone::Polycone{T})::T where T
    c = zero(T)
    for cone in pcone.sections
        c += capacity(cone)
    end
    return c
end

function surface(pcone::Polycone{T})::T where T
    N = length(pcone.sections)
    (; sections, Δϕ) = pcone
    annulus(rmin, rmax) = Δϕ * (rmax^2 - rmin^2) / 2
    surf = zero(T)
    for i in 1:N
        (; rmin1, rmax1, rmin2, rmax2, z) = sections[i]
        # conical (lateral) surfaces
        dmin, dmax = rmin2 - rmin1, rmax2 - rmax1
        surf += Δϕ * ((rmin1 + rmin2) / 2 * √(dmin^2 + 4z^2) + (rmax1 + rmax2) / 2 * √(dmax^2 + 4z^2))
        # phi faces (trapezoids of height 2z)
        Δϕ < T(2π) && (surf += 2z * ((rmax1 - rmin1) + (rmax2 - rmin2)))
    end
    # faces perpendicular to z: end caps and the steps between consecutive sections
    surf += annulus(sections[1].rmin1, sections[1].rmax1) + annulus(sections[N].rmin2, sections[N].rmax2)
    for i in 1:N-1
        lo, hi = sections[i], sections[i+1]
        overlap = annulus(max(lo.rmin2, hi.rmin1), max(min(lo.rmax2, hi.rmax1), max(lo.rmin2, hi.rmin1)))
        surf += annulus(lo.rmin2, lo.rmax2) + annulus(hi.rmin1, hi.rmax1) - 2overlap
    end
    return surf
end

function extent(pcone::Polycone{T})::Tuple{Point3{T},Point3{T}} where T
    N = length(pcone.sections)
    rmax::T = 0
    rmin::T = 0
    z::T = 0
    for i in eachindex(pcone.sections)
        cone = pcone.sections[i]
        cone.rmax1 > rmax && (rmax = cone.rmax1)
        cone.rmax2 > rmax && (rmax = cone.rmax2)
        cone.rmin1 < rmin && (rmin = cone.rmin1)
        cone.rmin2 < rmin && (rmin = cone.rmin2)
        z += 2 * cone.z
    end
    low, high = extent(Tube{T}(rmin, rmax, 0, pcone.ϕ₀, pcone.Δϕ))
    return (low + Vector3{T}(0, 0, pcone.zᵢ[1] - pcone.sections[1].z), high + Vector3{T}(0, 0, pcone.zᵢ[N] + pcone.sections[N].z))
end

function inside(pcone::Polycone{T}, point::Point3{T})::Int64  where T
    N = length(pcone.sections)
    x, y, z = point
    indexLow  = section_index(pcone, z - tolerance(T))
    indexHigh = section_index(pcone, z + tolerance(T))
    indexLow < 0 && indexHigh < 0 && return kOutside
    # beyond one end (the section at the other side can be thinner than the tolerance)
    indexLow < 0 && return inside(pcone.sections[indexHigh], point-Vector3{T}(0,0,pcone.zᵢ[indexHigh]))
    indexHigh < 0 && return inside(pcone.sections[indexLow], point-Vector3{T}(0,0,pcone.zᵢ[indexLow]))
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
Base.contains(pcone::Polycone{T}, p::Point3{T}) where T = inside(pcone, p) == kInside

function distance_to_in(pcone::Polycone{T}, point::Point3{T}, dir::Vector3{T})::T where T
    N = length(pcone.sections)
    z = point[3]
    dz = dir[3]
    distance = T(Inf)
    increment = dz > 0 ? 1 : -1
    abs(dz) < tolerance(T) && (increment = 0)
    index = section_index(pcone, z)
    index == -1 && (index = 1)
    index == -2 && (index = N)
    while true
        distance = distance_to_in(pcone.sections[index], point - Vector3{T}(0,0,pcone.zᵢ[index]), dir)
        (distance < T(Inf) || increment == 0) && break
        index += increment
        (index < 1 || index > N) && break
    end
    return distance
end

function safety_to_in(pcone::Polycone{T}, point::Point3{T})::T where T
    N = length(pcone.sections)
    (; sections, zᵢ) = pcone
    z = point[3]
    index = section_index(pcone, z)

    index == -1 && (index = 1)
    index == -2 && (index = N)

    minSafety = 0
    safety = safety_to_in(sections[index], point - Vector3{T}(0,0,zᵢ[index]))
    safety < tolerance(T) && return safety

    minSafety = safety
    zbase = zᵢ[index] + sections[index].z
    # going right
    for i in index+1:N
        dz = zᵢ[i] - sections[i].z - zbase
        dz >= minSafety && break
        safety = safety_to_in(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
        safety < minSafety && (minSafety = safety)
    end
    # going left
    if index > 1
        zbase = zᵢ[index-1] - sections[index-1].z
        for i in index-1:-1:1
            dz = zbase - zᵢ[i] + sections[i].z
            dz >= minSafety && break
            safety = safety_to_in(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
            safety < minSafety && (minSafety = safety)
        end
    end
    return minSafety
end

function normal(pcone::Polycone{T}, point::Point3{T}) where T
    N = length(pcone.sections)
    inside(pcone, point) != kSurface && return zero(Vector3{T})
    (; sections, zᵢ) = pcone
    x, y, z = point
    r = sqrt(x * x + y * y)
    tol = tolerance(T)/2
    n = zero(Vector3{T})
    count = 0
    # Lateral surfaces (conical and phi planes) of the sections touching the point
    for i in 1:N
        lp = point - Vector3{T}(0, 0, zᵢ[i])
        abs(lp[3]) <= sections[i].z + tol || continue
        inside(sections[i], lp) == kSurface || continue
        ln, lc = _lateral_normal(sections[i], lp)
        n += ln; count += lc
    end
    # Faces perpendicular to z: end planes and steps between sections. At a given z-plane the
    # point is on a face if its radius is covered by the section below or above but not by both.
    covers(rmin, rmax, strict) = strict ? (rmin + tol < r < rmax - tol) : (rmin - tol <= r <= rmax + tol)
    for i in 1:N
        sec = sections[i]
        ztop, zbot = zᵢ[i] + sec.z, zᵢ[i] - sec.z
        if abs(z - ztop) <= tol           # top face of section i
            covered_above = i < N && abs(zᵢ[i+1] - sections[i+1].z - ztop) <= tol &&
                            covers(sections[i+1].rmin1, sections[i+1].rmax1, true)
            if covers(sec.rmin2, sec.rmax2, false) && !covered_above
                n += Vector3{T}(0, 0, 1); count += 1
            end
        end
        if abs(z - zbot) <= tol           # bottom face of section i
            covered_below = i > 1 && abs(zᵢ[i-1] + sections[i-1].z - zbot) <= tol &&
                            covers(sections[i-1].rmin2, sections[i-1].rmax2, true)
            if covers(sec.rmin1, sec.rmax1, false) && !covered_below
                n += Vector3{T}(0, 0, -1); count += 1
            end
        end
    end
    count == 0 && return zero(Vector3{T})
    return normalize(n)
end

function safety_to_out(pcone::Polycone{T}, point::Point3{T})::T where T
    N = length(pcone.sections)
    inside(pcone, point) == kOutside && return -1
    (; sections, zᵢ) = pcone
    z = point[3]
    index = section_index(pcone, z)
    index < 0 && return -1
    safety = safety_to_out(sections[index], point - Vector3{T}(0,0,zᵢ[index]))
    minSafety = safety
    (minSafety == T(Inf) || minSafety < tolerance(T)) && return 0
    zbase = zᵢ[index] + sections[index].z
    # going right
    for i in index+1:N
        dz = zᵢ[i] - sections[i].z  - zbase
        dz >= minSafety && break
        safety = safety_to_in(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
        safety < minSafety && (minSafety = safety)
    end
    # going left
    if index > 1
        zbase = zbase = zᵢ[index-1] - sections[index-1].z
        for i in index-1:-1:1
            dz = zbase - zᵢ[i] + sections[i].z
            dz >= minSafety && break
            safety = safety_to_in(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
            safety < minSafety && (minSafety = safety)
        end
    end
    return minSafety
end

function distance_to_out(pcone::Polycone{T}, point::Point3{T}, dir::Vector3{T})::T where T
    N = length(pcone.sections)
    z = point[3]
    dz = dir[3]
    distance = T(Inf)
    if N == 1
        return distance_to_out(pcone.sections[1], point - Vector3{T}(0,0,pcone.zᵢ[1]), dir)
    end
    indexLow  = section_index(pcone, z - tolerance(T))
    indexHigh = section_index(pcone, z + tolerance(T))
    if indexLow < 0 && indexHigh < 0 
        return T(-1)
    elseif indexLow < 0 && indexHigh > 0
        index = indexHigh
        inside(pcone.sections[index], point - Vector3{T}(0,0,pcone.zᵢ[index])) == kOutside && return T(-1)
    elseif indexLow != indexHigh && indexLow > 0
        isin = inside(pcone.sections[indexLow], point - Vector3{T}(0,0,pcone.zᵢ[indexLow]))
        index = isin == kOutside ? indexHigh : indexLow
    else
        index = indexLow
    end

    if index < 0
        return zero(T)
    else
        inside(pcone.sections[index], point - Vector3{T}(0,0,pcone.zᵢ[index])) == kOutside && return T(-1)
    end

    distance = T(0)
    increment = dz > 0 ? 1 : -1
    abs(dz) < tolerance(T) && (increment = 0)
    pn = point
    istep = 0
    while true
        # section surface case
        if distance != 0 || istep < 2
            pn = point + distance * dir - Vector3{T}(0,0,pcone.zᵢ[index])
            inside(pcone.sections[index], pn ) == kOutside && break
        else
            pn -=  Vector3{T}(0,0,pcone.zᵢ[index])
        end
        istep += 1
        dist = distance_to_out(pcone.sections[index], pn, dir)
        dist == -1 && return distance
        if abs(dist) < tolerance(T)/2
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