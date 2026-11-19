# Geom4hep.jl

A native Julia geometry modeller for high energy physics, in the spirit of the
[Geant4](https://geant4.web.cern.ch) and [VecGeom](https://gitlab.cern.ch/VecGeom/VecGeom) geometry:
shapes (solids), boolean combinations, hierarchical volumes and materials, and a navigator that runs on
CPUs and GPUs.

Geom4hep.jl gives you:

- **Shapes** with the same parameters as the Geant4 solids (`Box`, `Tube`, `Cone`, `Trap`, `Sphere`,
  `Torus`, `Polycone`, `Polyhedra`, ...), boolean shapes and multi-unions, and the functions that a
  particle transport needs: [`inside`](@ref), [`distance_to_in`](@ref), [`distance_to_out`](@ref),
  [`safety_to_in`](@ref), [`safety_to_out`](@ref) and [`normal`](@ref).
- **Geometry models**: a hierarchy of logical volumes ([`Volume`](@ref)) placed inside each other, built
  in Julia or read from a [GDML](https://gdml.web.cern.ch) file with [`read_gdml`](@ref).
- **Navigation on CPUs and GPUs**: [`close_geometry`](@ref) turns a hierarchy of volumes into a flat
  [`GeometryModel`](@ref) of plain arrays, without dynamic dispatch, that a navigator ([`locate`](@ref),
  [`compute_step`](@ref)) traverses with the same code on the CPU and in GPU kernels written with
  [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) (Metal, CUDA, AMDGPU).
- **Visualization** of shapes and volumes with [Makie](https://docs.makie.org).
- **Benchmarks** of the shape functions with the methodology of VecGeom, to compare with it.

All the types are parametric on the floating point type, `Float64` or `Float32` (most GPUs, and Apple
GPUs in particular, are much faster, or only work, in single precision).

The geometry concepts are those of Geant4: see the chapter
[Geometry](https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geometry.html)
of the Geant4 Application Developer Guide for a description of the solids and of logical and physical
volumes.

## Installation

Geom4hep.jl is not yet registered. Install it from GitHub:

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/JuliaHEP/Geom4hep.jl")
```

The drawing functions are a package extension, activated when a Makie backend (`GLMakie`,
`CairoMakie` or `WGLMakie`) is loaded. To run on a GPU, load the KernelAbstractions backend of your GPU:
`Metal`, `CUDA` or `AMDGPU`.

## Quick start

This example builds a box of water with a lead ball inside, closes the geometry and follows a ray
through it. The [User Guide](@ref) explains each part.

```julia
using Geom4hep
const T = Float64

# Volumes: shape and material
world = Volume{T}("world", Box{T}(100, 100, 100), Material{T}("vacuum"; density=0))
water = Volume{T}("water", Box{T}(50, 50, 50), Material{T}("water"; density=1.0))
ball  = Volume{T}("ball", Sphere{T}(0, 20), Material{T}("lead"; density=11.35))

# Placements: the ball at (10, 0, 0) of the water box, the water box at the centre of the world
place_daughter!(water, Transformation3D{T}(10, 0, 0), ball)
place_daughter!(world, one(Transformation3D{T}), water)

# Flat model for navigation
model = close_geometry(world)

# Follow a ray from `point` along `dir` until it leaves the world
function trace(model, point, dir)
    state = locate(model, point)
    while in_world(state)
        step, next = compute_step(model, state, point, dir, T(Inf))
        println("step of ", step, " mm in material ", material_index(model, state))
        point += step * dir
        state = next
    end
end
trace(model, Point3{T}(-90, 0, 0), Vector3{T}(1, 0, 0))
```

## Status

The package is being revived and re-organized. The shapes most used in detector descriptions are
supported (see the [table of supported shapes](@ref Supported-shapes)); the full CMS 2018 GDML
geometry can be navigated, on CPUs and on Apple GPUs. Next steps are more solids (`Para`,
`ExtrudedSolid`, `TessellatedSolid`, ...), material properties and safety distances for boolean
shapes. See the [Release Notes](@ref).
