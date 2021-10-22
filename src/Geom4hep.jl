module Geom4hep

export Point3, Vector3, project
export AbstractShape, AbstractMaterial
export Box, Trd
export getindex, capacity, surface, extent, normal, distanceToOut, distanceToIn, inside, safetyToOut, safetyToIn, toMesh
export Material
export Transformation3D, RotMatrix3, RotXYZ, one, isone, transform, hasrotation, hastranslation, inv, lmul!
export Volume, placeDaughter!, draw
export NavigatorState, computeStep!, locateGlobalPoint!, currentVolume, getClosestDaughter
export kTolerance

using StaticArrays, GeometryBasics, LinearAlgebra, Rotations

# basic stuff
Vector3 = SVector{3}
Base.:*(x::Vector3, y::Vector3) = Vector3(x[1]*y[1], x[2]*y[2], x[3]*y[3])
Base.:-(p1::Point3, p2::Point3) = Vector3(p1[i]-p2[i] for i = 1:3)
project(p::Point3{T}, idx::Int64) where T<:AbstractFloat = Point3{T}([idx == i ? 0 : p[i] for i = 1:3]...)
const kTolerance = 1e-9

#= enums are neded to be exported
macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end
@exported_enum EInside 
=#
const kInside  = 0
const kSurface = 1
const kOutside = 2
export kInside, kSurface, kOutside

abstract type AbstractShape{T<:AbstractFloat} end
abstract type AbstractMaterial{T<:AbstractFloat} end

const coordmap = Dict(:dx => 1, :dy => 2, :dz => 3)

include("Transformation3D.jl")
include("Box.jl")
include("Trd.jl")
include("Volume.jl")
include("Navigators.jl")
include("Drawing.jl")


end # module
