#---Shapes---------------------------------------------------------------------------
const Shape{T<:AbstractFloat} = Union{BaseShape{T}, 
                                      BooleanUnion{T}, 
                                      BooleanSubtraction{T}, 
                                      BooleanIntersection{T}}
# const Shape{T<:AbstractFloat} = AbstractShape{T}
#---Volume-------------------------------------------------------------------------
struct VolumeP{T<:AbstractFloat,PV}
    id::Int64
    label::String
    shape::Shape{T}                     # Reference to the actual shape
    material::Material{T}               # Reference to material
    daughters::Vector{PV}
end

_last_volume_id = 1

function VolumeP{T, PV}(label::String, shape::Shape{T}, material::Material{T}, daughters::Vector{PV}) where {T<:AbstractFloat, PV}
    global _last_volume_id
    _last_volume_id += 1
    VolumeP(_last_volume_id, label, shape, material, daughters)
end


#---PlacedVolume-------------------------------------------------------------------
struct PlacedVolume{T<:AbstractFloat}
    idx::Int64
    transformation::Transformation3D{T}
    volume::VolumeP{T,PlacedVolume{T}}
    aabb::AABB{T}
end
PlacedVolume(idx::Int64,transformation::Transformation3D{T},volume::VolumeP{T,PlacedVolume{T}}) where {T} = PlacedVolume{T}(idx,transformation,volume)
function PlacedVolume{T}(idx::Int64,transformation::Transformation3D{T},volume::VolumeP{T,PlacedVolume{T}}) where {T}
    if volume.shape isa NoShape
        aabb=AABB{T}(Point3{T}(0,0,0),Point3{T}(0,0,0))
    else
        aabb=AABB{T}(transform_extent(extent(volume.shape),transformation)...)
    end
    return PlacedVolume{T}(idx,transformation,volume,aabb)
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

@inline function distanceToIn(pvol::PlacedVolume{T}, p::Point3{T}, d::Vector3{T}, rcp_d::Vector3{T}=inv.(d))::T where T<:AbstractFloat
    !intersect(pvol.aabb, p, d, rcp_d) && return T(Inf)
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
        box = Box{T}((high - low)/2. * 1.05)  # increase the bounding box
        tra = Transformation3D{T}(one(RotMatrix3{T}), -(high + low)/2.)
        mat = Material{T}("vacuum"; density=0)
        world = Volume{T}("world", box, mat)
        placeDaughter!(world, tra, vol)
        return world
    end
end


function AABB(pvol::PlacedVolume{T}) where T
    AABB{T}(transform_extent(extent(pvol.volume.shape),pvol.transformation)...)
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
    tm=typemax(T)
    mmin = Point3{T}(tm,tm,tm)

    tm=typemin(T)
    mmax = Point3{T}(tm,tm,tm)
    for pvol in agg.pvolumes
        pmin = pvol.aabb.min
        pmax = pvol.aabb.max
        mmin = min.(mmin,pmin)
        mmax = min.(max,pmax)
    end

    return (mmin,mmax)
end

function inside(agg::Aggregate{T}, point::Point3{T})::Int64  where T<:AbstractFloat
    for pvol in agg.pvolumes
        inout = inside(pvol.volume.shape, pvol.transformation * point)
        (inout == kInside || inout == kSurface) && return inout
    end
    return kOutside
end

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
