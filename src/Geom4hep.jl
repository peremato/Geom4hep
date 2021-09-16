module Geom4hep

export Point, Vec
export Box, getindex, capacity, surface, extent, normal

using StaticArrays, GeometryBasics, LinearAlgebra

abstract type AbstractShape{T<:AbstractFloat} end

const coordmap = Dict(:dx => 1, :dy => 2, :dz => 3)

include("Box.jl")

end # module
