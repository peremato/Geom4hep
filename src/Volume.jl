#---Used for Aggregates------------------------------------------------------------- 
struct NoShape{T,PV} <: AbstractShape{T}
    pvolumes::Vector{PV}
end
NoShape{T}() where T = NoShape{T,Nothing}([])

#---Shape--------------------------------------------------------------------------
const Shape{T} = Union{NoShape{T},
                       Box{T},
                       Trd{T},
                       Trap{T},
                       Tube{T},
                       Cone{T},
                       Polycone{T},
                       CutTube{T},
                       BooleanUnion{T},
                       BooleanIntersection{T}, 
                       BooleanSubtraction{T}} where T<:AbstractFloat

#---Volume-------------------------------------------------------------------------
struct VolumeP{T<:AbstractFloat,PV}
    label::String
    shape::Shape{T}                     # Reference to the actual shape
    material::Material{T}               # Reference to material
    daughters::Vector{PV}
end

#---PlacedVolume-------------------------------------------------------------------
struct PlacedVolume{T<:AbstractFloat}
    idx::Int64
    transformation::Transformation3D{T}
    volume::VolumeP{T,PlacedVolume{T}}
end

#---Convenient Alias to simplify signatures---------------------------------------
const Volume{T} = VolumeP{T,PlacedVolume{T}} where T<:AbstractFloat  

#---Constructor--------------------------------------------------------------------
function Volume{T}(label::String, shape::Shape{T}, material::Material{T}) where T<:AbstractFloat
    Volume{T}(label, shape, material, Vector{PlacedVolume{T}}())   # call the default constructor
end

#---Utilities----------------------------------------------------------------------
function Base.show(io::IO, vol::Volume{T}) where T
    name = vol.label
    print(io, "Volume{$T} name = $name")
end

function AbstractTrees.children(vol::Volume{T}) where T
    collect((d.volume for d in vol.daughters))
end

function Base.getindex(vol::Volume{T}, indx...) where T
    v = vol
    for i in indx
        v = v.daughters[i].volume
    end
    return v
end

function contains(pvol::PlacedVolume{T}, p::Point3{T})::Bool where T<:AbstractFloat
    inside(pvol.volume.shape, pvol.transformation * p) == kInside
end
function contains(vol::Volume{T}, p::Point3{T})::Bool where T<:AbstractFloat
    inside(vol.shape, p) == kInside
end
function distanceToIn(pvol::PlacedVolume{T}, p::Point3{T}, d::Vector3{T})::T where T<:AbstractFloat
    distanceToIn(pvol.volume.shape, pvol.transformation * p, pvol.transformation * d)
end


function placeDaughter!(volume::Volume{T}, placement::Transformation3D{T}, subvol::Volume{T}) where T<:AbstractFloat
    if subvol.shape isa NoShape
        for d in subvol.daughters
            push!(volume.daughters, PlacedVolume(length(volume.daughters)+1, d.transformation * placement, d.volume))
        end
    else
        push!(volume.daughters, PlacedVolume(length(volume.daughters)+1, placement,subvol))
    end
    volume
end

function getWorld(vol::Volume{T}) where T<:AbstractFloat
    if vol.label == "world"
        return vol
    else
        low, high = extent(vol.shape)
        box = Box{T}((high - low)/2. .+ 0.1)  # increase the bounding box
        tra = Transformation3D{T}(one(RotMatrix3{T}), -(high + low)/2.)
        mat = Material{T}("vacuum"; density=0)
        world = Volume{T}("world", box, mat)
        placeDaughter!(world, tra, vol)
        return world
    end
end

#---Bounding Boxes---------------------------------------------------------------------------------
# Axis Aligned Bounding Box, defined by the min/max coordinates
struct AABB{T<:AbstractFloat}
    min::Point3{T}
    max::Point3{T}
end

_area(d) = 2(d[1]*d[2] + d[2]*d[3] + d[3]*d[1])
_volume(d) = d[1]*d[2]*d[3]

diagonal(bb::AABB{T}) where T = bb.max .- bb.min
area(bb::AABB{T}) where T = _area(diagonal(bb))
volume(bb::AABB{T}) where T = _volume(diagonal(bb))

function GeometryBasics.coordinates(bb::AABB{T}) where T
    a, b = bb.min, bb.max
    Point3{T}[a, (b[1], a[2], a[3]), (a[1], b[2], a[3]), (a[1], a[2], b[3]),
                 (a[1], b[2], b[3]), (b[1], a[2], b[3]), (b[1], b[2], a[3]), b]
end
GeometryBasics.faces(::AABB{T}) where T = QuadFace{Int64}[(1,3,7,2), (1,2,6,4), (1,4,5,3), (8,6,2,7), (8,7,3,5), (8,5,4,6)]
GeometryBasics.normals(::AABB{T}) where T = Point3{T}[(0,0,-1), (0,-1,0), (-1,0,0), (1,0,0), (0,1,0), (0,0,1)]

inside(bb::AABB{T}, point::Point3{T}) where T =  all(bb.min .< point .< bb.max)
function intersect(bb::AABB{T}, point::Point3{T},dir::Vector3{T}) where T
    point = point - (bb.max + bb.min)/2
    dimens = (bb.max - bb.min)/2
    distsurf = Inf
    distance = -Inf
    distout = Inf
    for i in 1:3
        din  = (-copysign(dimens[i],dir[i]) - point[i])/dir[i]
        tout =   copysign(dimens[i],dir[i]) - point[i]
        dout = tout/dir[i]
        dsur = copysign(tout, dir[i])
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
    (distance >= distout || distout <= kTolerance(T)/2 || abs(distsurf) <= kTolerance(T)/2) ? false : true
end

function AABB(pvol::PlacedVolume{T}) where T
    a, b = extent(pvol.volume.shape)
    min_ = [ Inf,  Inf,  Inf]
    max_ = [ -Inf,  -Inf,  -Inf] 
    coords =  Point3{T}[a, (b[1], a[2], a[3]), (a[1], b[2], a[3]), (a[1], a[2], b[3]),
                           (a[1], b[2], b[3]), (b[1], a[2], b[3]), (b[1], b[2], a[3]), b]
    for c in coords
        p = c * pvol.transformation
        for i in 1:3 
            p[i] > max_[i] && (max_[i] = p[i])
            p[i] < min_[i] && (min_[i] = p[i])
        end
    end
    return AABB{T}(Point3{T}(min_), Point3{T}(max_))
end


#---Assemblies--------------------------------------------------------------------------------------
struct Assembly{T<:AbstractFloat}
end

function Assembly{T}(label::String) where T<:AbstractFloat
    pvolumes = Vector{PlacedVolume{T}}()
    Volume{T}(label, NoShape{T,PlacedVolume{T}}(pvolumes), Material{T}("vacuum"; density=0), pvolumes)   # call the default constructor
end

const Aggregate{T} =  NoShape{T,PlacedVolume{T}} where T

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, shape::NoShape{T,PlacedVolume{T}}) where T
    print(io, "Aggregate{$T} with ", length(shape.pvolumes), " volumes: \n")
    for d in shape.pvolumes
        print(io, "  $d\n")
    end
end

function GeometryBasics.coordinates(agg::Aggregate{T}, facets=36) where {T<:AbstractFloat}
    return (point * pvol.transformation for pvol in agg.pvolumes for point in coordinates(pvol.volume.shape, facets))
end

function GeometryBasics.faces(agg::Aggregate{T}, facets=36) where {T<:AbstractFloat}
    offset = [length(coordinates(pvol.volume.shape, facets)) for pvol in agg.pvolumes]
    accumulate!(+, offset, offset; init=-offset[1])
    indexes = Vector{TriangleFace{Int}}()
    for (i, pvol) in enumerate(agg.pvolumes)
        append!(indexes, [t .+ offset[i] for t in faces(pvol.volume.shape, facets)])
    end
    return indexes
end

#---Basic Functions--------------------------------------------------------------------------------
function capacity(agg::Aggregate{T}) where T<:AbstractFloat
    sum((capacity(pvol.volume.shape) for pvol in agg.pvolumes))
end

function surface(agg::Aggregate{T}) where T<:AbstractFloat
    sum((surface(pvol.volume.shape) for pvol in agg.pvolumes))
end

function extent(agg::Aggregate{T})::Tuple{Point3{T},Point3{T}} where T<:AbstractFloat
    return (min((extent(pvol.volume.shape)[1] * pvol.transformation  for pvol in agg.pvolumes)...),
            max((extent(pvol.volume.shape)[2] * pvol.transformation  for pvol in agg.pvolumes)...))
end

# FIXME never used
# function inside(agg::Aggregate{T}, point::Point3{T})::Int64  where T<:AbstractFloat
#     for pvol in agg.pvolumes
#         inout = inside(pvol.volume.shape, pvol.transformation * point)
#         (inout == kInside || inout == kSurface) && return inout
#     end
#     return kOutside
# end

function safetyToOut(agg::Aggregate{T}, point::Point3{T}) where T<:AbstractFloat
end

function safetyToIn(agg::Aggregate{T}, point::Point3{T}) where T<:AbstractFloat
end

function distanceToOut(agg::Aggregate{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    for pvol in agg.pvolumes
        lpoint = pvol.transformation * point
        inout = inside(pvol.volume.shape, lpoint)
        inout == kSurface && return T(0)
        inout == kInside && return distanceToOut(pvol.volume.shape, lpoint, pvol.transformation * dir)
    end
    return T(-1)
end

function distanceToIn(agg::Aggregate{T}, point::Point3{T}, dir::Vector3{T})::T where T<:AbstractFloat
    distance = T(Inf)
    for i in eachindex(agg.pvolumes)
        pvol = agg.pvolumes[i]
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
