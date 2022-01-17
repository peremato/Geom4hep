#---Used for Agregates 
struct NoShape{T} <: AbstractShape{T}
    shapename::String
end

#---Shape--------------------------------------------------------------------------
const Shape{T} = Union{NoShape{T},Box{T},Trd{T},Tube{T},Cone{T}} where T<:AbstractFloat

#---Volume-------------------------------------------------------------------------
struct Mother{T,PV}
    label::String
    shape::Shape{T}                     # Reference to the actual shape
    material::Material                  # Reference to material
    daughters::Vector{PV} 
end

#---PlacedVolume-------------------------------------------------------------------
struct PlacedVolume{T<:AbstractFloat}
    transformation::Transformation3D{T}
    volume::Mother{T,PlacedVolume{T}}
end

const Volume{T} = Mother{T,PlacedVolume{T}} where T<:AbstractFloat

contains(pvol::PlacedVolume{T}, p::Point3{T}) where T<:AbstractFloat = inside(pvol.volume.shape, pvol.transformation * p) == kInside
contains(vol::Volume{T}, p::Point3{T}) where T<:AbstractFloat = inside(vol.shape, p) == kInside
distanceToIn(pvol::PlacedVolume{T}, p::Point3{T}, d::Vector3{T}) where T<:AbstractFloat = distanceToIn(pvol.volume.shape, pvol.transformation * p, pvol.transformation * d)

function Volume{T}(label::String, shape::Shape{T}, material::Material) where T<:AbstractFloat
    Volume{T}(label, shape, material, Vector{PlacedVolume{T}}())
end

function placeDaughter!(volume::Volume{T}, placement::Transformation3D{T}, subvol::Volume{T}) where T<:AbstractFloat
    if subvol.shape isa NoShape
        for d in subvol.daughters
            push!(volume.daughters, PlacedVolume(placement * d.transformation, d.volume))
        end
    else
        push!(volume.daughters, PlacedVolume(placement,subvol))
    end
end


