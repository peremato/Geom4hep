module Geom4hep

export Point3, Vector3, DPoint3, move!
export AbstractShape, AbstractMaterial
export Box, Trd, TBox, TTrd
export getindex, capacity, surface, extent, normal, distanceToOut,  distanceToIn, inside, safetyToOut, safetyToIn, toMesh
export Material, Isotope, Element
export Transformation3D, RotMatrix3, RotXYZ, one, isone, transform, hasrotation, hastranslation, inv, lmul!
export Mother, PlacedVolume, Volume, placeDaughter!, draw
export NavigatorState, computeStep!, locateGlobalPoint!, reset!, currentVolume, getClosestDaughter
export kTolerance
export processGDML
export Triangle, Intersection, intersect, distanceToPlane, PV

using StaticArrays, GeometryBasics, LinearAlgebra, Rotations

include("BasicTypes.jl")
include("Transformation3D.jl")
include("TriangleIntersect.jl")
include("Trd.jl")
include("Box.jl")
include("Materials.jl")
include("Volume.jl")
include("Navigators.jl")
include("Drawing.jl")
include("GDML.jl")

end # module
