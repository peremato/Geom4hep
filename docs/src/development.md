# Development

## Structure of the package

| Directory | Content |
|:----------|:--------|
| `src/shapes` | One file per shape: its type, the shape functions and the mesh for drawing |
| `src/model` | The flat [`GeometryModel`](@ref): its types, the CSG programs of boolean shapes and [`close_geometry`](@ref) |
| `src/navigation` | The navigator ([`locate`](@ref), [`compute_step`](@ref)) and [`xray`](@ref) |
| `src/GDML.jl` | The GDML reader |
| `src/Benchmark.jl` | The shape benchmarks |
| `ext/` | The Makie extension (drawing) |
| `test/` | The unit tests |
| `examples/` | Examples with their own environment, and GDML files |
| `docs/` | This documentation |

## Running the tests

```bash
julia --project=. -e 'import Pkg; Pkg.test()'
```

## Building the documentation

```bash
julia --project=docs -e 'import Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The examples are converted from their Literate sources (`docs/src/examples/*.lit`) and executed by
`docs/make.jl`; the generated pages, notebooks and scripts are not committed. Open
`docs/build/index.html` to see the result. The documentation is built and deployed by the CI for each
push to master and for each pull request (as a preview).

## Adding a shape

1. Add the shape type and its functions (`inside`, `distance_to_in`, `distance_to_out`,
   `safety_to_in`, `safety_to_out`, `normal`, `extent`, `capacity`, `surface`, and `coordinates` and
   `faces` for the mesh) in a new file in `src/shapes`, and include it in `src/Geom4hep.jl`.
2. Add it to `BaseShape` in `src/shapes/Boolean.jl`, and to the flat model in
   `src/model/GeometryModel.jl` and `src/model/CloseGeometry.jl`. Shapes with variable-size data (as
   polycones) keep it in shared arrays of the model.
3. Read it from GDML in `src/GDML.jl`.
4. Add tests, comparing with Geant4 or VecGeom where possible, also in `Float32`.
5. Check the size of the GPU kernels: Apple's compiler inlines all calls, and fails (or takes very
   long) on too large kernels. Avoid adding call sites of the shape functions and unrolled loops
   that call them (see the `@nounroll` macro).
