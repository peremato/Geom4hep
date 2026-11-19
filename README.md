# Geom4hep
A native Julia geometry modeller for HEP, in the spirit of the Geant4/VecGeom geometry: basic solids, boolean combinations, hierarchical volumes and materials, with data structures designed to be usable on GPUs.

Current functionality:

- Basic shapes: `Box`, `Tube`, `Trd`, `Trap`, `Cone`, `CutTube`, `Sphere` (also GDML `orb`), `Torus`, `Polycone`, `Polyhedra`
- Boolean shapes: `BooleanUnion`/`BooleanSubtraction`/`BooleanIntersection`, and `MultiUnion` (union of any number of placed shapes)
- Materials: `Isotope`, `Element`, `Material`
- Shape functions: `inside`, `distance_to_in`, `distance_to_out`, `safety_to_in`, `safety_to_out`, `normal`, `extent`, `capacity`, `surface`
- Construction of hierarchical geometry models with `Volume`, `PlacedVolume` and `Assembly`
- GDML reader to input geometries (`read_gdml`)
- System of units in `Geom4hep.Units` (`using Geom4hep.Units; Box{Float64}(1cm, 2cm, 3cm)`)
- Flat geometry model for navigation on CPUs and GPUs (see below)
- 3D graphics support for shapes and geometry models with [Makie](https://docs.makie.org), loaded automatically as a package extension when a Makie backend is loaded
- Benchmarks of the shape functions following the methodology of VecGeom, and comparison with VecGeom

All types are parametric on the floating point type (`Float64` or `Float32`).

## Supported shapes
The solids of Geant4 and their status in Geom4hep (✅ supported, ❌ not yet supported):

| Geant4 solid | Geom4hep | GDML element | Notes |
|---|---|---|---|
| `G4Box` | ✅ `Box` | `box` | |
| `G4Tubs` | ✅ `Tube` | `tube` | |
| `G4CutTubs` | ✅ `CutTube` | `cutTube` | |
| `G4Cons` | ✅ `Cone` | `cone` | |
| `G4Trd` | ✅ `Trd` | `trd` | |
| `G4Trap` | ✅ `Trap` | `trap` | |
| `G4Sphere` | ✅ `Sphere` | `sphere` | with ϕ and θ segments |
| `G4Orb` | ✅ `Sphere` | `orb` | as a full `Sphere` (no dedicated, faster type) |
| `G4Torus` | ✅ `Torus` | `torus` | |
| `G4Polycone` | ✅ `Polycone` | `polycone` | z-planes form only; not the (r, z) corners form (`genericPolycone`) |
| `G4Polyhedra` | ✅ `Polyhedra` | `polyhedra` | z-planes form only; not the (r, z) corners form (`genericPolyhedra`) |
| `G4UnionSolid`, `G4SubtractionSolid`, `G4IntersectionSolid` | ✅ `BooleanUnion`, `BooleanSubtraction`, `BooleanIntersection` | `union`, `subtraction`, `intersection` | safeties are conservative (0) in the flat model |
| `G4MultiUnion` | ✅ `MultiUnion` | ❌ `multiUnion` | no acceleration structure over the components |
| `G4Para` | ❌ | `para` | |
| `G4GenericTrap` | ❌ | `arb8` | |
| `G4GenericPolycone` | ❌ | `genericPolycone` | |
| `G4EllipticalTube` | ❌ | `eltube` | |
| `G4EllipticalCone` | ❌ | `elcone` | |
| `G4Ellipsoid` | ❌ | `ellipsoid` | |
| `G4Paraboloid` | ❌ | `paraboloid` | |
| `G4Hype` | ❌ | `hype` | |
| `G4Tet` | ❌ | `tet` | |
| `G4ExtrudedSolid` | ❌ | `xtru` | |
| `G4TessellatedSolid` | ❌ | `tessellated` | |
| `G4TwistedBox`, `G4TwistedTrd`, `G4TwistedTrap`, `G4TwistedTubs` | ❌ | `twistedbox`, `twistedtrd`, `twistedtrap`, `twistedtubs` | |
| `G4ScaledSolid` | ❌ | `scaledSolid` | |
| `G4ReflectedSolid` | ❌ | `reflectedSolid` | |

GDML solids that are not supported are replaced by `NoShape` with a warning.

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
shared by CPU and GPU. Next: more solids (para, extruded, tessellated, ...), material
properties, and safety distances for boolean shapes.

Known limitations:
- Other GDML solids (`para`, `xtru`, `tessellated`, ...) are not yet supported; they are replaced by `NoShape` with a warning.
- Polyhedra: only the z-planes form (GDML `polyhedra`), not the (r, z) corners form (`genericPolyhedra`).
- In Float32, the tolerance (1 µm) is coarse: a ray running parallel to a surface, within the tolerance,
  can be taken as inside or outside the shape (in the X-ray images of the CMS geometry, a few pixels in
  10⁵ differ from Float64).
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

### Benchmarks of the shapes
`benchmark_shape` times the shape functions (`inside`, `contains`, `distance_to_in`, `safety_to_in`,
`distance_to_out`, `safety_to_out`) with the methodology of the shape benchmarks of VecGeom (world box,
sampling of the points and directions, 50% of the points inside, 80% of the rays hitting the shape),
and `vecgeom_suite()` has the setups of the VecGeom benchmark programs (`test/benchmark/*Benchmark.cpp`)
for the shapes of Geom4hep. The example runs the suite and reports the time per call, directly on the
shapes and through the dispatch of the flat model:
```
julia --project=examples examples/Benchmarks.jl [name filter]
```
With `--vecgeom DIR` it also runs the VecGeom benchmark programs found in `DIR` with the same numbers of
points and repetitions, and compares the timings (VecGeom configured with `-DBUILD_TESTING=ON
-DVECGEOM_TEST_BENCHMARK=ON`; the programs are in `<build>/test/benchmark`). The results are also
written to a CSV file (`--csv FILE`, default `benchmarks.csv`).

### Draw a Detector
The example `DrawDetector.jl` constructs a detector model and draws a 3D image of it. To keep the graphics window alive run it with the julia option `-i` (interactive).
```
julia --project=examples -i examples/DrawDetector.jl
```

### Produce a X-Ray image of a detector
The example `XRay.jl` generates X-Ray images of a detector model (part of the CMS 2018 geometry): the mass thickness traversed along straight lines parallel to the axes, computed on the CPU or on a GPU (`Metal` is part of the examples environment; add `CUDA` to it for NVIDIA GPUs). The result is saved in `cms2018-xray-<device>.png`.
```
julia --project=examples -t auto examples/XRay.jl [cpu|metal|cuda]
```

### Density slices of a detector
The example `Render.jl` draws the density of the material in slices of the CMS 2018 geometry (`cms2018-slices.png`).
```
julia --project=examples -t auto examples/Render.jl
```
