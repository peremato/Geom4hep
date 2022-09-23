module Geom4hep

export Point3, Point2, Vector3, Vector2, nonzero
export AbstractShape, AbstractMaterial
export Shape, NoShape, Box, Trd, TBox, TTrd, Tube, Wedge, isInside, isOutside, Cone, Polycone, getNz, getNSections, getSectionIndex
export CutTube, Plane, Boolean, Trap, BooleanUnion
export getindex, capacity, surface, extent, normal, distanceToOut,  distanceToIn, inside, safetyToOut, safetyToIn
export Material, Isotope, Element
export Transformation3D, RotMatrix3, RotXYZ, one, isone, transform, hasrotation, hastranslation, inv, lmul!
export PlacedVolume, Volume, Assembly, placeDaughter!, draw, draw!, drawDistanceToIn, drawDistanceToOut, getWorld, children
export AABB, area, volume, BVHParam, buildBVH, BVH, pvolindices, pushPvolIndices!
export AbstractNavigator, TrivialNavigator, BVHNavigator, NavigatorState, computeStep!, locateGlobalPoint!, reset!, isInVolume, currentVolume, getClosestDaughter, containedDaughters
export kTolerance
export processGDML
export Triangle, Intersection, intersect, distanceToPlane
export Tesselation, coordinates, faces, normals, mesh

using Requires
using StaticArrays, GeometryBasics, LinearAlgebra, Rotations, AbstractTrees
include("BasicTypes.jl")
include("Transformation3D.jl")
include("TriangleIntersect.jl")
include("Trd.jl")
include("Box.jl")
include("Wedge.jl")
include("Tube.jl")
include("Cone.jl")
include("Polycone.jl")
include("CutTube.jl")
using Unityper
include("Boolean.jl")
include("Trap.jl")
include("Materials.jl")
include("Volume.jl")
include("BVH.jl")
include("Navigators.jl")
include("GDML.jl")
include("./DistanceIn.jl")
include("./DistanceOut.jl")
include("./Distance.jl")
include("./Inside.jl")
include("Benchmark.jl")
#include("CuGeom.jl")  # PackageCompiler fails?

function __init__()
    @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("Drawing.jl")
end

end # module
