# User Guide

This guide describes how to build a geometry with Geom4hep.jl and how to navigate it. The concepts
are those of the Geant4 geometry (solids, logical volumes and placements); the
[Geometry](https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geometry.html)
chapter of the Geant4 Application Developer Guide describes them in more detail. The
[Examples](@ref Running-the-examples) show all of this running, with their output.

## Floating point type and units

All the types are parametric on the floating point type `T`, `Float64` or `Float32`, which is given
explicitly when creating objects: `Box{Float32}(1, 2, 3)`. The parameters are converted to `T`. A
geometry must use a single type; [`read_gdml`](@ref) takes it as an argument.

`Float32` is needed on most GPUs (Apple GPUs do not support `Float64`). With `Float32`, the surfaces
are thicker ([`tolerance`](@ref) is 1 µm instead of 1 nm), so that rounding errors do not let a ray
cross a surface without seeing it.

Lengths are in mm and angles in radians. The module `Geom4hep.Units` defines units with the Geant4
conventions:

```julia
using Geom4hep, Geom4hep.Units
box  = Box{Float64}(1cm, 2cm, 3cm)          # half lengths 10, 20 and 30 mm
tube = Tube{Float64}(0, 5cm, 1m, 0, 90deg)  # a quarter of a cylinder
```

## Points, vectors and transformations

[`Point3`](@ref) is a position and [`Vector3`](@ref) a displacement or a direction. Both are static
vectors (from StaticArrays.jl), so they support indexing, broadcasting, `dot`, `cross`, `norm`, ...
The difference of two points is a vector, and a point plus a vector is a point. Directions given to
the shape functions must be unit vectors.

```julia
p = Point3{Float64}(1, 2, 3)
d = normalize(Vector3{Float64}(1, 1, 0))     # using LinearAlgebra
q = p + 10 * d                               # a Point3
```

A [`Transformation3D`](@ref) places a daughter volume in its mother (or the right operand of a boolean
shape in the frame of the left one). Applied to a point of the mother frame, it gives the point in
the frame of the daughter: `t * p` is `R * (p - t)` for a translation `t` and rotation `R`.
`Transformation3D{T}(x, y, z)` therefore places the daughter at (x, y, z). Rotations come from
[Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl):

```julia
t1 = Transformation3D{Float64}(10, 0, 0)                          # translation
t2 = Transformation3D{Float64}(0, 0, 50, RotXYZ{Float64}(0, 0, π/4)) # rotation and translation
t3 = Transformation3D{Float64}(0, 0, 0, π/2, 0, 0)               # RotXYZ angles
p * t1      # inverse transformation: from the daughter frame to the mother frame
t1 * t2     # composition
```

Note the convention, the same as Geant4 and VecGeom use internally: the rotation `R` rotates points of
the mother frame into the daughter frame, so the daughter appears rotated by the inverse of `R`.

## Shapes

A shape (solid) is defined in its own frame, usually centred at the origin. The parameters are the
same as those of the corresponding Geant4 solid (half lengths, radii, and ϕ and θ ranges with a
start angle and an extent).

### Supported shapes

The Geant4 solids and their status in Geom4hep (✅ supported, ❌ not yet supported):

| Geant4 solid | Geom4hep | GDML element | Notes |
|:-------------|:---------|:-------------|:------|
| `G4Box` | ✅ [`Box`](@ref) | `box` | |
| `G4Tubs` | ✅ [`Tube`](@ref) | `tube` | |
| `G4CutTubs` | ✅ [`CutTube`](@ref) | `cutTube` | |
| `G4Cons` | ✅ [`Cone`](@ref) | `cone` | |
| `G4Trd` | ✅ [`Trd`](@ref) | `trd` | |
| `G4Trap` | ✅ [`Trap`](@ref) | `trap` | |
| `G4Sphere` | ✅ [`Sphere`](@ref) | `sphere` | with ϕ and θ segments |
| `G4Orb` | ✅ [`Sphere`](@ref) | `orb` | as a full `Sphere` |
| `G4Torus` | ✅ [`Torus`](@ref) | `torus` | |
| `G4Polycone` | ✅ [`Polycone`](@ref) | `polycone` | z-planes form only |
| `G4Polyhedra` | ✅ [`Polyhedra`](@ref) | `polyhedra` | z-planes form only |
| `G4UnionSolid`, `G4SubtractionSolid`, `G4IntersectionSolid` | ✅ [`BooleanUnion`](@ref), [`BooleanSubtraction`](@ref), [`BooleanIntersection`](@ref) | `union`, `subtraction`, `intersection` | safeties are conservative (0) in the flat model |
| `G4MultiUnion` | ✅ [`MultiUnion`](@ref) | ❌ `multiUnion` | no acceleration structure over the components |
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
| `G4TwistedBox`, `G4TwistedTrd`, `G4TwistedTrap`, `G4TwistedTubs` | ❌ | `twistedbox`, ... | |
| `G4ScaledSolid` | ❌ | `scaledSolid` | |
| `G4ReflectedSolid` | ❌ | `reflectedSolid` | |

The constructors:

```julia
T = Float64
Box{T}(x, y, z)                                   # half lengths
Tube{T}(rmin, rmax, z, ϕ₀, Δϕ)                    # z is the half length
Cone{T}(rmin1, rmax1, rmin2, rmax2, z, ϕ₀, Δϕ)    # radii at -z and +z
CutTube{T}(rmin, rmax, z, ϕ₀, Δϕ, bottom_normal, top_normal)
Trd{T}(x1, x2, y1, y2, z)                         # half lengths at -z and +z
Trap{T}(z, θ, ϕ, y1, x1, x2, α1, y2, x3, x4, α2)  # as G4Trap
Sphere{T}(rmin, rmax, ϕ₀, Δϕ, θ₀, Δθ)
Torus{T}(rmin, rmax, rtor, ϕ₀, Δϕ)
Polycone{T}(rmin, rmax, z, ϕ₀, Δϕ)                # vectors of radii at the z planes
Polyhedra{T}(nsides, rmin, rmax, z, ϕ₀, Δϕ)       # radii to the sides (apothems)
```

### Shape functions

Every shape implements the functions that a particle transport needs. The points and directions are
in the frame of the shape:

| Function | Result |
|:---------|:-------|
| [`inside(shape, p)`](@ref inside) | `kInside`, `kSurface` or `kOutside` |
| [`distance_to_in(shape, p, d)`](@ref distance_to_in) | distance along `d` from an outside point to the shape (`Inf` if it misses it) |
| [`distance_to_out(shape, p, d)`](@ref distance_to_out) | distance along `d` from an inside point to the surface |
| [`safety_to_in(shape, p)`](@ref safety_to_in) | lower bound of the distance from an outside point to the shape |
| [`safety_to_out(shape, p)`](@ref safety_to_out) | lower bound of the distance from an inside point to the surface |
| [`normal(shape, p)`](@ref normal) | outward unit normal at a point of the surface |
| [`extent(shape)`](@ref extent) | `(lower, upper)` corners of a bounding box |
| [`capacity(shape)`](@ref capacity), [`surface(shape)`](@ref surface) | volume and surface area |

The surfaces have a thickness: a point closer than `tolerance(T)/2` to a surface is on it. A point on
the surface moving into the shape gets `distance_to_in == 0`, and moving out of it
`distance_to_out == 0`, so that a navigator that stops on a surface does not get stuck there. The
safeties are cheaper than the exact distances: they are used to take steps in all directions (e.g.
for multiple scattering) without computing intersections.

### Boolean shapes

[`BooleanUnion`](@ref), [`BooleanSubtraction`](@ref) and [`BooleanIntersection`](@ref) combine two
shapes; the right one is placed in the frame of the left one. They can be nested:

```julia
box  = Box{T}(10, 10, 10)
hole = Tube{T}(0, 5, 20, 0, 2π)
s1 = BooleanSubtraction(box, hole)                                        # a box with a hole along z
s2 = BooleanUnion(s1, Sphere{T}(0, 5), Transformation3D{T}(0, 0, 10))     # and a dome on top
```

A [`MultiUnion`](@ref) is the union of any number of placed shapes:

```julia
mu = MultiUnion{T}()
for i in 1:10
    add_component!(mu, Box{T}(1, 1, 1), Transformation3D{T}(3i, 0, 0))
end
```

## Materials

A [`Material`](@ref) has a name, a state, a density and a composition in [`Element`](@ref)s, which
are made of [`Isotope`](@ref)s. Only the density is used for now (by [`xray`](@ref)):

```julia
water = Material{T}("water"; density=1.0)
```

## Volumes and placements

A [`Volume`](@ref) (Geant4's logical volume) is a shape filled with a material. Daughter volumes are
placed in it with [`place_daughter!`](@ref) (Geant4's physical volume placement), which records a
[`PlacedVolume`](@ref) with the transformation. The same volume can be placed many times, in the same
or in different mothers. The daughters must be inside the mother and must not overlap each other;
this is not checked.

```julia
world  = Volume{T}("world", Box{T}(1000, 1000, 1000), Material{T}("air"; density=0.0012))
layer = Volume{T}("module", Box{T}(100, 100, 10), Material{T}("iron"; density=7.87))
for i in -2:2
    place_daughter!(world, Transformation3D{T}(0, 0, 100i), layer)
end
world[3]          # the volume of the 3rd daughter
```

An [`Assembly`](@ref) groups volumes without a mother shape: placing an assembly places each of its
volumes. [`world_volume`](@ref) returns a world volume that contains a given volume.

## GDML

[`read_gdml(file, T)`](@ref read_gdml) reads a [GDML](https://gdml.web.cern.ch) file, the geometry
exchange format of Geant4 and ROOT, and returns its world volume:

```julia
world = read_gdml("cms2018.gdml", Float32)
```

The materials, the shapes in the table of [supported shapes](@ref Supported-shapes), the booleans,
the volumes, the physical volumes and the assemblies are read, with the units of the file converted to
mm and radians. Solids that are not supported are replaced by [`NoShape`](@ref) with a warning.

## Geometry model and navigation

### Closing the geometry

The hierarchy of volumes is convenient to build a geometry, but not to navigate it: its objects have
different types, are referenced by pointers and cannot be moved to a GPU.
[`close_geometry`](@ref) converts it to a [`GeometryModel`](@ref):

- all the data is in arrays of plain (isbits) structures, referenced by indices: the model can be
  moved to a GPU with `Adapt.adapt`;
- the shapes are stored by type, and a shape function is called with [`shape_call`](@ref) through a
  generated `if` chain over the types instead of dynamic dispatch;
- boolean shapes, of any nesting, become flat CSG programs over primitive shapes, and long chains of
  unions become multi-unions;
- the daughters of volumes with many daughters are organized in a bounding volume hierarchy (BVH).

```julia
model = close_geometry(world)
model, index = close_geometry(world; with_index=true)   # index.volumes[i] is the Volume of volume i
```

### Navigating

A [`NavState`](@ref) is the position of a point in the hierarchy: the path of placements from the
world to the current volume, and the transformation from the global frame to the frame of the current
volume. It is immutable, so it can be kept in registers on a GPU.

- [`locate(model, point)`](@ref locate) gives the state of a global point (the deepest volume that
  contains it).
- [`compute_step(model, state, point, dir, step_limit)`](@ref compute_step) gives the distance to the
  next boundary along `dir` (or `step_limit`) and the state after crossing it.
- [`in_world(state)`](@ref in_world) is `false` once the point has left the world, and
  [`material_index(model, state)`](@ref material_index) gives the material of the current volume.

```julia
function track_length(model, point, dir)
    state = locate(model, point)
    length = zero(eltype(point))
    while in_world(state)
        step, state = compute_step(model, state, point, dir, typemax(eltype(point)))
        point += step * dir
        length += step
    end
    return length
end
```

After crossing a boundary, the new state is computed for the point pushed by
[`push_tolerance`](@ref) beyond it. Like Geant4, the navigator remembers the daughter that a step has
just left ("blocked" volume), to avoid re-entering it immediately because of rounding errors.

### GPUs

The navigator functions are plain Julia functions on isbits data, so they can be called from GPU
kernels. [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) gives kernels
that run on any GPU. [`xray`](@ref) is an example: it computes an X-ray image of a geometry (the mass
thickness along parallel rays) on the CPU (with Julia threads) or on a GPU:

```julia
world = read_gdml("cms2018.gdml", Float32)
model = close_geometry(world)
u, v, image, nsteps = xray(model, 3, 1000)                                  # CPU
using Metal
u, v, image, nsteps = xray(model, 3, 1000; backend=MetalBackend())          # Apple GPU
```

Some GPU compilers (Apple's in particular) inline all the function calls of a kernel: the size of a
kernel grows with the number of shape types of the model and with the number of call sites of the
shape functions. This is why the model dispatches without recursion and the CSG programs are flat.

## Visualization

The drawing functions need a Makie backend (`GLMakie` for interactive windows, `CairoMakie` for
images). [`draw`](@ref) creates a figure with a shape or a volume hierarchy, and [`draw!`](@ref) adds it
to an existing `LScene`:

```julia
using GLMakie: GLMakie, Figure, LScene
draw(Torus{Float64}(1, 2, 10, 0, 3π/2))
fig = Figure(); scene = LScene(fig[1, 1])
draw!(scene, world; wireframe=true, maxlevel=3)
```

Makie exports some names that are also Geom4hep shapes (`Box`, `Sphere`, ...): import only the Makie
names you need, as above, or qualify the shapes (`Geom4hep.Box`).

[`mesh(shape)`](@ref mesh) returns a `GeometryBasics.Mesh` of the surface of a shape, for other
visualization or export tools.

## Benchmarks

[`benchmark_shape`](@ref) times the shape functions of a shape with the methodology of the VecGeom
shape benchmarks: the shape in a world box, points sampled in the box (half of them inside the shape)
and rays of which 80% hit the shape. [`vecgeom_suite`](@ref) has the setups of the VecGeom benchmark
programs for the shapes of Geom4hep. The script `examples/Benchmarks.jl` runs the suite and, given a
VecGeom build, runs the VecGeom programs with the same setups and compares the times:

```
julia --project=examples examples/Benchmarks.jl --vecgeom <vecgeom-build>/test/benchmark Polycone
```
