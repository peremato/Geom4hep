#---Material-----------------------------------------------------------------------
abstract type AbstractShape{T<:AbstractFloat} end
abstract type AbstractMaterial{T<:AbstractFloat} end

struct Material{T<:AbstractFloat} <: AbstractMaterial{T}
    label::String
    density::T
end
Material(label::String, density::Float64) = Material{Float64}(label,density)


#---Volume-------------------------------------------------------------------------
abstract type AbstractPlacedVolume{T<:AbstractFloat} end
struct Volume{T<:AbstractFloat}
    label::String
    shape::Geom4hep.AbstractShape{T}    # Reference to the actual shape
    material::AbstractMaterial{T}       # Reference to material
    daughters::Vector{AbstractPlacedVolume{T}}
end

#---PlacedVolume-------------------------------------------------------------------
struct PlacedVolume{T<:AbstractFloat} <: AbstractPlacedVolume{T}
    transformation::Transformation3D{T}
    volume::Volume{T}
end

contains(pvol::PlacedVolume{T}, p::Point3{T}) where T<:AbstractFloat = inside(pvol.volume.shape, pvol.transformation * p) == kInside
contains(vol::Volume{T}, p::Point3{T}) where T<:AbstractFloat = inside(vol.shape, p) == kInside
distanceToIn(pvol::PlacedVolume{T}, p::Point3{T}, d::Vector3{T}) where T<:AbstractFloat = distanceToIn(pvol.volume.shape, pvol.transformation * p, pvol.transformation * d)

function Volume(label::String, shape::Geom4hep.AbstractShape{T}, material::AbstractMaterial{T}) where T<:AbstractFloat
    Volume{T}(label, shape, material, [])
end

function placeDaughter!(volume::Volume{T}, placement::Transformation3D{T}, subvol::Volume{T}) where T<:AbstractFloat
    push!(volume.daughters, PlacedVolume(placement,subvol))
end


