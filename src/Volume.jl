#---Used for Agregates 
struct NoShape{T} <: AbstractShape{T}
end

#---Shape--------------------------------------------------------------------------
const Shape{T} = Union{NoShape{T},
                       Box{T},
                       Trd{T},
                       Trap{T},
                       Tube{T},
                       Cone{T},
                       Polycone{T},
                       CutTube{T},
                       Boolean{T}} where T<:AbstractFloat

#---Volume-------------------------------------------------------------------------
struct VolumeP{T,PV}
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

contains(pvol::PlacedVolume{T}, p::Point3{T}) where T<:AbstractFloat = inside(pvol.volume.shape, pvol.transformation * p) == kInside
contains(vol::Volume{T}, p::Point3{T}) where T<:AbstractFloat = inside(vol.shape, p) == kInside
distanceToIn(pvol::PlacedVolume{T}, p::Point3{T}, d::Vector3{T}) where T<:AbstractFloat = distanceToIn(pvol.volume.shape, pvol.transformation * p, pvol.transformation * d)


function placeDaughter!(volume::Volume{T}, placement::Transformation3D{T}, subvol::Volume{T}) where T<:AbstractFloat
    if subvol.shape isa NoShape
        for d in subvol.daughters
            push!(volume.daughters, PlacedVolume(length(volume.daughters)+1, placement * d.transformation, d.volume))
        end
    else
        push!(volume.daughters, PlacedVolume(length(volume.daughters)+1, placement,subvol))
    end
end

function getWorld(vol::Volume{T}) where T<:AbstractFloat
    if vol.label == "world"
        return vol
    else
        low, high = extent(vol.shape)
        box = Box{T}((high - low)/2. .+ kTolerance(T)*10)
        tra = Transformation3D{T}(one(RotMatrix3{T}), -(high + low)/2.)
        mat = Material{T}("vacuum"; density=0)
        world = Volume{T}("world", box, mat)
        placeDaughter!(world, tra, vol)
        return world
    end
end
