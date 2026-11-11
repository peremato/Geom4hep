module Geom4hep

export Point3, Point2, Vector3, Vector2, nonzero
export AbstractShape, AbstractMaterial
export Shape, NoShape, Box, Trd, Tube, Wedge, isInside, isOutside, Cone, Polycone, getNz, getNSections, getSectionIndex
export CutTube, Plane, Trap, BooleanUnion, BooleanIntersection, BooleanSubtraction
export getindex, capacity, surface, extent, normal, distanceToOut,  distanceToIn, inside, safetyToOut, safetyToIn
export Material, Isotope, Element
export Transformation3D, RotMatrix3, RotXYZ, one, isone, transform, hasrotation, hastranslation, inv, lmul!
export PlacedVolume, Volume, Assembly, placeDaughter!, draw, draw!, drawDistanceToIn, drawDistanceToOut, getWorld, children
export AABB, area, volume, BVHParam, buildBVH, BVH, pvolindices
export AbstractNavigator, TrivialNavigator, BVHNavigator, NavigatorState, computeStep!, locateGlobalPoint!, reset!, isInVolume, currentVolume, getClosestDaughter, containedDaughters
export kTolerance, kPushTolerance
export processGDML
export Tessellation, coordinates, faces, normals, mesh

using StaticArrays, GeometryBasics, LinearAlgebra, Rotations, AbstractTrees, DataStructures
using GeometryBasics: mesh, Mesh
include("Units.jl")
include("BasicTypes.jl")
include("Transformation3D.jl")
include("BoundingBox.jl")
#---Shapes (order matters: helpers first, Boolean last)
include("shapes/Wedge.jl")
include("shapes/Plane.jl")
include("shapes/Trd.jl")
include("shapes/Box.jl")
include("shapes/Tube.jl")
include("shapes/Cone.jl")
include("shapes/Polycone.jl")
include("shapes/CutTube.jl")
include("shapes/Trap.jl")
include("shapes/NoShape.jl")
include("shapes/Boolean.jl")
include("Materials.jl")
include("Volume.jl")
include("BVH.jl")
include("Navigators.jl")
include("GDML.jl")
include("Benchmark.jl")
#include("CuGeom.jl")  # CUDA-only prototype, to be replaced by a KernelAbstractions implementation

#---Functions implemented in package extensions (see ext/)------------------------------------------
"""
    draw(object; kwargs...)

Draw a shape, volume hierarchy, AABB or BVH. Requires a Makie backend (e.g. `using GLMakie`).
"""
function draw end
"""
    draw!(scene, object; kwargs...)

Draw a shape, volume hierarchy, AABB or BVH into an existing `LScene`. Requires a Makie backend.
"""
function draw! end
function drawDistanceToIn end
function drawDistanceToOut end

end # module
