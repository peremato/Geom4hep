# Geom4hep
A native Julia geometry modeller for HEP, in the spirit of the Geant4/VecGeom geometry: basic solids, boolean combinations, hierarchical volumes and materials, with data structures designed to be usable on GPUs.

Current functionality:

- Basic shapes: `Box`, `Tube`, `Trd`, `Trap`, `Cone`, `CutTube`, `Polycone`
- Boolean shapes: `BooleanUnion`/`BooleanSubtraction`/`BooleanIntersection`, and `MultiUnion` (union of any number of placed shapes)
- Materials: `Isotope`, `Element`, `Material`
- Shape functions: `inside`, `distance_to_in`, `distance_to_out`, `safety_to_in`, `safety_to_out`, `normal`, `extent`, `capacity`, `surface`
- Construction of hierarchical geometry models with `Volume`, `PlacedVolume` and `Assembly`
- GDML reader to input geometries (`read_gdml`)
- System of units in `Geom4hep.Units` (`using Geom4hep.Units; Box{Float64}(1cm, 2cm, 3cm)`)
- Flat geometry model for navigation on CPUs and GPUs (see below)
- 3D graphics support for shapes and geometry models with [Makie](https://docs.makie.org), loaded automatically as a package extension when a Makie backend is loaded
- Benchmark functions for shapes and a number of examples

All types are parametric on the floating point type (`Float64` or `Float32`).

## Navigation on CPUs and GPUs
`close_geometry(world)` converts a hierarchy of volumes into a flat `GeometryModel`: all the data is stored
in arrays of isbits structures referenced by indices, shapes are called without dynamic dispatch, boolean
shapes become flat CSG programs over primitive shapes, and the daughters of each volume are organized in a
bounding volume hierarchy. The navigator works on this model with immutable states, so the same code runs
on the CPU and in GPU kernels ([KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl)):

```julia
using Geom4hep
world = read_gdml("examples/cms2018.gdml", Float32)
model = close_geometry(world)
point, dir = Point3{Float32}(0, 0, 0), Vector3{Float32}(1, 0, 0)
state = locate(model, point)                                         # deepest volume containing the point
step, state = compute_step(model, state, point, dir, Float32(Inf))   # distance to the next boundary

# X-ray image (mass thickness along rays parallel to the z axis), on the CPU or on a GPU
u, v, image, nsteps = xray(model, 3, 1000)                          # CPU (Julia threads)
using Metal; u, v, image, nsteps = xray(model, 3, 1000; backend=MetalBackend())   # Apple GPU
```

The older tree-based navigators (`TrivialNavigator`, `BVHNavigator`) are still available but have known
problems (they can miss boundaries of non-convex mother volumes) and will be removed.

## Status
The package is being revived and re-organized. Done: package hygiene and tests, own point/vector types and
a `snake_case` API, consistent normals and Float32 support, and the flat geometry model with a navigator
shared by CPU and GPU. Next: more solids (polyhedra, sphere, torus, tessellated, ...), material
properties, and safety distances for boolean shapes.

Known limitations:
- `polyhedra`, `torus` and other GDML solids are not yet supported; they are replaced by `NoShape` with a warning.
- In Float32, the tolerance (1 µm) is coarse for very thin features: a small fraction of the distances to
  boolean shapes with sub-tolerance details can be wrong (about 0.03% for the CMS geometry).
- The CPU backend of KernelAbstractions 0.10 (POCL) is not used: `xray` runs on Julia threads on the CPU.

## Running the unit tests
```
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Examples
The examples have their own environment in `examples/Project.toml` (with GLMakie, BenchmarkTools, ...). Instantiate it once from the repository root:
```
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```
(on Julia 1.10 also run `julia --project=examples -e 'using Pkg; Pkg.develop(path=".")'`).

### Benchmarks
To run the predefined set of benchmarks for shapes do:
```
julia --project=examples examples/Benchmarks.jl
```
Note that the times are the total time for all points being considered. Divide by the number of points to obtain the execution time per function invocation.

### Draw a Detector
The example `DrawDetector.jl` constructs a detector model and draws a 3D image of it. To keep the graphics window alive run it with the julia option `-i` (interactive).
```
julia --project=examples -i examples/DrawDetector.jl
```

### Produce a X-Ray image of a detector
The example `XRay.jl` generates X-Ray images of a detector model (part of the CMS 2018 geometry): the mass thickness traversed along straight lines parallel to the axes, computed on the CPU or on a GPU (add `Metal` or `CUDA` to the examples environment first). The result is saved in `cms2018-xray.png`.
```
julia --project=examples -t auto examples/XRay.jl [cpu|metal|cuda]
```

### Density slices of a detector
The example `Render.jl` draws the density of the material in slices of the CMS 2018 geometry (`cms2018-slices.png`).
```
julia --project=examples -t auto examples/Render.jl
```
