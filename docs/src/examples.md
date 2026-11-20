# Running the examples

The pages in this section are generated from the [Literate](https://fredrikekre.github.io/Literate.jl)
sources in [`docs/src/examples`](https://github.com/JuliaHEP/Geom4hep.jl/tree/master/docs/src/examples)
and run with the current version of Geom4hep.jl when the documentation is built. They show the code
together with its output and plots. Each page has links to download the example as a Jupyter notebook
or a plain Julia script.

| Example | Content |
|:--------|:--------|
| [Shapes](examples/Shapes.md) | Creating shapes, calling the shape functions and drawing them |
| [Volumes](examples/Volumes.md) | Building a hierarchy of volumes and drawing it |
| [Navigation](examples/Navigation.md) | Closing a geometry, following rays through it and X-ray images |
| [GDML](examples/GDML.md) | Reading GDML files: the CMS 2018 detector |
| [RoundCube](examples/RoundCube.md) | Porting a Geant4 solid built from boolean operations, and checking it |
| [Benchmarks](examples/Benchmarks.md) | Timing the shape functions with the methodology of VecGeom |

## Setting up an environment

The examples use extra packages for plotting (`CairoMakie` or `GLMakie`). The simplest way to run them
is with the environment of this documentation, which lists every package they need:

```bash
git clone https://github.com/JuliaHEP/Geom4hep.jl.git
cd Geom4hep.jl
julia --project=docs -e 'import Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/src/examples/Shapes.jl
```

To run an example downloaded from these pages instead, put it in a directory with the
[`Project.toml`](https://github.com/JuliaHEP/Geom4hep.jl/blob/master/docs/Project.toml) of the
documentation, run `julia --project=. -e 'import Pkg; Pkg.add(url="https://github.com/JuliaHEP/Geom4hep.jl"); Pkg.instantiate()'`,
and then the example, or open the notebook in Jupyter.

The repository has more examples in the [`examples`](https://github.com/JuliaHEP/Geom4hep.jl/tree/master/examples)
directory, with their own environment: X-ray images on GPUs (`XRay.jl`), density slices of the CMS
detector (`Render.jl`), interactive drawing of a detector (`DrawDetector.jl`) and the comparison of the
shape benchmarks with VecGeom (`Benchmarks.jl`).
