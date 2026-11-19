# Public API

Documentation of the public interface of `Geom4hep.jl`.

## Index
```@index
Pages = ["api.md"]
```

## Points, vectors and transformations
```@docs
Point3
Vector3
Point2
Transformation3D
```

## Shapes
```@docs
AbstractShape
Box
Trd
Trap
Tube
Cone
CutTube
Sphere
Torus
Polycone
Polyhedra
BooleanUnion
BooleanSubtraction
BooleanIntersection
MultiUnion
add_component!
NoShape
```

## Shape functions
```@docs
inside
distance_to_in
distance_to_out
safety_to_in
safety_to_out
normal
extent
capacity
surface
estimate_capacity
estimate_surface
tolerance
push_tolerance
```

## Materials
```@docs
AbstractMaterial
Isotope
Element
Material
```

## Volumes
```@docs
Volume
PlacedVolume
place_daughter!
Assembly
world_volume
read_gdml
```

## Geometry model and navigation
```@docs
close_geometry
GeometryModel
ShapeRef
shape_call
NavState
locate
compute_step
in_world
material_index
xray
```

## Visualization
```@docs
draw
draw!
draw_distance_to_in
draw_distance_to_out
mesh
```

## Benchmarks
```@docs
ShapeBenchmark
benchmark_shape
vecgeom_suite
```

## Units
```@docs
Geom4hep.Units
```
