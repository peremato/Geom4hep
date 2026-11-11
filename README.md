# Geom4hep
A native Julia geometry modeller for HEP, in the spirit of the Geant4/VecGeom geometry: basic solids, boolean combinations, hierarchical volumes and materials, with data structures designed to be usable on GPUs.

Current functionality:

- Basic shapes: `Box`, `Tube`, `Trd`, `Trap`, `Cone`, `CutTube`, `Polycone`, and `BooleanUnion`/`BooleanSubtraction`/`BooleanIntersection`
- Materials: `Isotope`, `Element`, `Material`
- Essential navigation functions on shapes (`distance_to_in`, `distance_to_out`, `safety_to_in`, `safety_to_out`, `inside`, ...)
- Construction of hierarchical geometry models with `Volume`, `PlacedVolume` and `Assembly`
- Very basic navigator `TrivialNavigator` (linear search among all daughters)
- Bounding Volume Hierarchy accelerated navigator `BVHNavigator`
- GDML reader to input geometries (`read_gdml`)
- System of units in `Geom4hep.Units` (`using Geom4hep.Units; Box{Float64}(1cm, 2cm, 3cm)`)
- 3D graphics support for shapes and geometry models with [Makie](https://docs.makie.org), loaded automatically as a package extension when a Makie backend is loaded
- Benchmark functions for shapes and a number of examples

All types are parametric on the floating point type (`Float64` or `Float32`).

## Status
The package is being revived and re-organized. Planned steps: package hygiene and tests (in progress), source re-layout with own point/vector types and a `snake_case` API, a flat index-based geometry model shared by CPU and GPU navigation (via KernelAbstractions.jl), and more solids and material properties.

Known limitations:
- `polyhedra`, `torus` and other GDML solids are not yet supported; they are replaced by `NoShape` with a warning.
- The CUDA prototype (`src/CuGeom.jl` and the `CUDA*` examples) is not loaded and will be replaced by a KernelAbstractions implementation.

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
The example `XRay.jl` generates X-Ray images of a detector model (part of the CMS 2018 geometry). It uses the density of the traversed material along straight lines parallel to the axes, and saves the result in `cms2018.png`.
```
julia --project=examples examples/XRay.jl
```
