using Revise
using Geom4hep
using StaticArrays
using GeometryBasics

#---Cone--------------------------------------------------------------------------------------------
struct Polycone{T<:AbstractFloat, N, L} <: AbstractShape{T}
    sections::SVector{N, Cone{T}}  # array of sections
    zᵢ::SVector{N,T}               # sections shifts
    rmin::SVector{L,T}             # rmin Vector
    rmax::SVector{L,T}             # rmax Vector
    z::SVector{L,T}                # z-planes
    ϕ₀::T                          # starting ϕ value (in radians)
    Δϕ::T                          # delta ϕ value of tube segment (in radians)
end

#---Constructor------------------------------------------------------------------------------------
function Polycone{T}(rmin::Vector, rmax::Vector, z::Vector, ϕ₀, Δϕ) where T<:AbstractFloat
    @assert length(rmin) == length(rmax) == length(z) "Vectors for rmin, rmax and z are different size"
    L = length(rmin)
    sections = [Cone{T}(rmin[i],rmax[i], rmin[i+1], rmax[i+1], (z[i+1]-z[i])/2, ϕ₀, Δϕ) for i in 1:L-1 if (z[i+1] - z[i]) > kTolerance(T)]
    zoffs = [(z[i]+z[i+1])/2 for i in 1:L-1 if (z[i+1] - z[i]) > kTolerance(T)]
    N = length(sections)
    Polycone{T,N,L}(sections, zoffs, rmin, rmax, z, ϕ₀, Δϕ)
end

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, pcone::Polycone{T,N}) where {T,N}
    print(io, "Polycone{$T, $N}",(rmin=pcone.rmin, rmax=pcone.rmax, z=pcone.z, ϕ₀=pcone.ϕ₀, Δϕ=pcone.Δϕ))
end

function GeometryBasics.coordinates(pcone::Polycone{T,N}, facets=36) where {T,N}
    return (coord + Vector3{T}(0,0,pcone.zᵢ[i]) for i in 1:N for coord in coordinates(pcone.sections[i], facets))
end

function GeometryBasics.faces(pcone::Polycone{T,N}, facets=36) where {T,N}
    (; sections, Δϕ) = pcone
    infacets = facets
    issector = Δϕ < 2π
    issector ?  facets =  round(Int64, (facets/2π) * Δϕ) : nothing
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    nbf = Int(facets / 2)            # Number of faces
    nbv = issector ? nbf + 1 : nbf   # Number of vertices
    offset = fill(0,N)
    for i in 1:N-1
        ishollow = sections[i].rmin1 > 0 || sections[i].rmin2 > 0
        nbc = ishollow ? nbv : 1         # Number of centers
        offset[i+1] = offset[i] + 2 * nbv + 2 * nbc
    end
    @show offset, N
    indexes = Vector{TriangleFace{Int}}()
    for i in 1:N
        append!(indexes, [t .+ offset[i] for t in faces(pcone.sections[i], infacets)])
        @show length(faces(pcone.sections[i]))
    end
    return indexes
end

#---Basic functionality-----------------------------------------------------------------------------
getNz(::Polycone{T,N}) where {T,N} = N + 1
getNSections(::Polycone{T,N}) where {T,N} = N
function getSectionIndex(pcone::Polycone{T,N}, zpos) where {T,N}
    (; sections, zᵢ) = pcone
    zpos < zᵢ[1] - sections[1].z && return -1
    for i  in 1:N
        Δz = sections[i].z
        z  = zᵢ[i]
        zpos >= z - Δz && zpos <= z + Δz && return i
    end
    return -2
end

function capacity(pcone::Polycone{T,N})::T where {T,N}
    sum(map(capacity, pcone.sections))
end

function surface(pcone::Polycone{T,N})::T where {T,N}
    surf = T(0)
    for i in 1:N
        (; rmin1, rmax1, rmin2, rmax2, z, Δϕ) = pcone.sections[i]  
        mmin = (rmin1 + rmin2) * 0.5
        mmax = (rmax1 + rmax2) * 0.5
        dmin = (rmin2 - rmin1)
        dmax = (rmax2 - rmax1)
        surf += Δϕ * (mmin * √(dmin * dmin + 4 * z * z) + mmax * √(dmax * dmax + 4 * z * z))
        if i == 1
            surf +=  Δϕ * 0.5 * (rmax1 * rmax1 - rmin1 * rmin1)
        elseif i == N
            surf +=  Δϕ * 0.5 * (rmax2 * rmax2 - rmin2 * rmin2)
        end
    end
    return surf
end

function extent(pcone::Polycone{T,N})::Tuple{Point3{T},Point3{T}} where {T,N}
    rmax::T = 0
    rmin::T = 0
    z::T = 0
    for cone in pcone.sections
        cone.rmax1 > rmax && (rmax = cone.rmax1)
        cone.rmax2 > rmax && (rmax = cone.rmax2)
        cone.rmin1 < rmin && (rmin = cone.rmin1)
        cone.rmin2 < rmin && (rmin = cone.rmin2)
        z += 2 * cone.z
    end
    low, high = extent(Tube{T}(rmin, rmax, 0, pcone.ϕ₀, pcone.Δϕ))
    return (low + Vector3{T}(0, 0, pcone.zᵢ[1] - pcone.sections[1].z), high + Vector3{T}(0, 0, pcone.zᵢ[N] + pcone.sections[N].z))
end

function inside(pcone::Polycone{T,N}, point::Point3{T}) where {T,N}
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
Base.contains(pcone::Polycone{T,N}, p::Point3{T}) where {T,N} = inside(pcone, p) == kInside

function distanceToIn(pcone::Polycone{T,N}, point::Point3{T}, dir::Vector3{T})::T where {T,N}
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

function safetyToIn(pcone::Polycone{T,N}, point::Point3{T})::T where {T,N}
    (; sections, zᵢ) = pcone
    z = point[3]
    index = getSectionIndex(pcone, z)

    index == -1 && (index = 1)
    index == -2 && (index = N)

    minSafety = 0
    safety = safetyToIn(sections[index], point - Vector3{T}(0,0,zᵢ[index]))
    safety < kTolerance(T) && return safety

    minSafety = safety
    zbase = zᵢ[index] + sections[index].z
    # going right
    for i in index+1:N
        dz = zᵢ[i] - sections[i].z - zbase
        dz >= minSafety && break
        safety = safetyToIn(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
        safety < minSafety && (minSafety = safety)
    end
    # going left
    if index > 1
        zbase = zᵢ[index-1] - sections[index-1].z
        for i in index-1:-1:1
            dz = zbase - zᵢ[i] + sections[i].z
            dz >= minSafety && break
            safety = safetyToIn(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
            safety < minSafety && (minSafety = safety)
        end
    end
    return minSafety
end

function safetyToOut(pcone::Polycone{T,N}, point::Point3{T})::T where {T,N}
    inside(pcone, point) == kOutside && return -1
    (; sections, zᵢ) = pcone
    z = point[3]
    index = getSectionIndex(pcone, z)
    index < 0 && return -1
    safety = safetyToOut(sections[index], point - Vector3{T}(0,0,zᵢ[index]))
    minSafety = safety
    (minSafety == Inf || minSafety < kTolerance(T)) && return 0
    zbase = zᵢ[index] + sections[index].z
    # going right
    for i in index+1:N
        dz = zᵢ[i] - sections[i].z  - zbase
        dz >= minSafety && break
        safety = safetyToIn(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
        safety < minSafety && (minSafety = safety)
    end
    # going left
    if index > 1
        zbase = zbase = zᵢ[index-1] - sections[index-1].z
        for i in index-1:-1:1
            dz = zbase - zᵢ[i] + sections[i].z
            dz >= minSafety && break
            safety = safetyToIn(sections[i], point - Vector3{T}(0,0,zᵢ[i]))
            safety < minSafety && (minSafety = safety)
        end
    end
    return minSafety
end

function distanceToOut(pcone::Polycone{T,N}, point::Point3{T}, dir::Vector3{T})::T where {T,N}
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
        # section surface case
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