module Geom4hep

export Point3, Point2, Vector3, Vector2, nonzero
export AbstractShape, AbstractMaterial
export Shape, NoShape, Box, Trd, TBox, TTrd, Tube, Wedge, isInside, isOutside, Cone
export getindex, capacity, surface, extent, normal, distanceToOut,  distanceToIn, inside, safetyToOut, safetyToIn
export Material, Isotope, Element
export Transformation3D, RotMatrix3, RotXYZ, one, isone, transform, hasrotation, hastranslation, inv, lmul!
export Mother, PlacedVolume, Volume, placeDaughter!, draw, draw!, drawDistanceToIn, drawDistanceToOut
export NavigatorState, computeStep!, locateGlobalPoint!, reset!, currentVolume, getClosestDaughter
export kTolerance
export processGDML
export Triangle, Intersection, intersect, distanceToPlane, PV
export Tesselation, coordinates, faces, normals

using StaticArrays, GeometryBasics, LinearAlgebra, Rotations

include("BasicTypes.jl")
include("Transformation3D.jl")
include("TriangleIntersect.jl")
include("Trd.jl")
include("Box.jl")
include("Wedge.jl")
include("Tube.jl")
include("Cone.jl")
include("Materials.jl")
include("Volume.jl")
include("Navigators.jl")
include("Drawing.jl")
include("GDML.jl")
include("Benchmark.jl")
include("CuGeom.jl")

end # module
