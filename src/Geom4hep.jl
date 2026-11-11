module Geom4hep

export Point3, Point2, Vector3, Vector2, nonzero
export AbstractShape, AbstractMaterial
export Shape, NoShape, Box, Trd, Tube, Wedge, is_inside, is_outside, Cone, Polycone, nz, nsections, section_index
export CutTube, Plane, Trap, BooleanUnion, BooleanIntersection, BooleanSubtraction
export getindex, capacity, surface, extent, normal, distance_to_out,  distance_to_in, inside, safety_to_out, safety_to_in
export Material, Isotope, Element
export Transformation3D, RotMatrix3, RotXYZ, one, isone, transform, hasrotation, hastranslation, inv, lmul!
export PlacedVolume, Volume, Assembly, place_daughter!, draw, draw!, draw_distance_to_in, draw_distance_to_out, world_volume, children
export AABB, area, volume, BVHParam, build_bvh, BVH, pvolindices
export AbstractNavigator, TrivialNavigator, BVHNavigator, NavigatorState, compute_step!, locate_global_point!, reset!, is_in_volume, current_volume, closest_daughter, contained_daughters
export tolerance, push_tolerance
export read_gdml
export coordinates, faces, normals, mesh

using StaticArrays, LinearAlgebra, Rotations, AbstractTrees, DataStructures
using GeometryBasics: GeometryBasics, Mesh, mesh, coordinates, faces, normals, TriangleFace, QuadFace
include("Units.jl")
include("Vectors.jl")
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
function draw_distance_to_in end
function draw_distance_to_out end

end # module
