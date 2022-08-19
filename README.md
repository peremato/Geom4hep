# Geom4hep
An exercise for evaluation of the Julia language for a Geometry modeling package used in HEP (VecGeom). The current functionality:

- Basic shapes: `Box`, `Tube`, `Trd`, `Cone`, `CutTube`, `Polycone`, `Boolean` 
- Materials: `Isotope`, `Element`, `Material`
- Essential functions to navigate on them (`distanceToIn`, `distanceToOut`, ...)
- Construction of hierarchical geometry models with `Volume` and `PlacedVolumes` 
- Very basic navigator `TrivialNavigator` (linear search among all daughters)
- Bounding Volume Hierarchy accelerated navigator `BVHNavigator`
- 3D graphics support for displaying the shapes and geometry models using the [Makie](https://makie.juliaplots.org/stable/) package
- GDML reader to input geometries
- Benchmark functions for shapes
- A number of examples to demonstrate the available functionality 
- Tests with using GPUs (CUDA.jl) 

## Running the unit tests
A number of unit tests is provided under /test for each of the geometrical shapes. To run all of them do:

```
julia --project=. test/runtests.jl 
```

## Examples

### Benchmarks
To run the predefined set of benchmarks for shapes do:
```
julia --project=. examples/Benchmarks.jl 
```
Note that the times are the total time for all points being consirered. Divide by the number of points to obtain the execution time per function invocation.

### Draw a Detector
The example `DrawDetector.jl` constructs detector model and draws a 3D image of it. In oder to keep the graphics window alive to be able to interact with it you need to run the example with the julia option `-i` (interactive).
```
julia --project=. -i  examples/DrawDetector.jl 
```

### Produce a X-Ray image of a detector
The example `XRay.jl` generates X-Ray images starting from a detector model (e.g. `world`). It defines a function that generates the data using the density of the traversed material from a detector model. You need to run the example with the julia option `-i` (interactive).
```
julia --project=. -i  examples/XRay.jl 
```

