#---Material-----------------------------------------------------------------------
abstract type AbstractShape{T<:AbstractFloat} end
abstract type AbstractMaterial{T<:AbstractFloat} end

struct Material{T<:AbstractFloat} <: AbstractMaterial{T}
    label::String
end
Material(label::String) = Material{Float64}(label)



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
    tranformation::Transformation3D{T}
    volume::Volume{T}
end


function Volume(l::String, s::Geom4hep.AbstractShape{T}, m::AbstractMaterial{T} ) where T <: AbstractFloat
    Volume{T}(l, s, m,[])
end
function placeDaughter!(volume::Volume{T}, placement::Transformation3D{T}, subvol::Volume{T}) where T <: AbstractFloat
    push!(volume.daughters, PlacedVolume(placement,subvol))
end

using GLMakie
function draw(vol::Volume) 
    scene = mesh(toMesh(vol.shape), color=:yellow, transparency=true, ambient=0.7, visible= vol.label == "world" ? false : true)
    for daughter in vol.daughters
        m = toMesh(daughter.volume.shape)
        points = GeometryBasics.coordinates(m)
        faces  = GeometryBasics.faces(m)
        map!(c -> daughter.tranformation * c, points, points)
        mesh!(points, faces, color=:blue, transparency=true, ambient=0.7)
    end
    display(scene)
end
