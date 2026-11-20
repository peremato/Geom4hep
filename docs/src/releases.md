# Release Notes

## 0.1.0 - 19-11-2026
- New features:
    - Shapes: `Sphere` (with ϕ and θ segments, also GDML `orb`), `Polyhedra` and `Torus`, and the mesh
      of all the shapes for drawing.
    - Flat geometry model (`close_geometry`) and navigator (`locate`, `compute_step`) shared by CPUs
      and GPUs (KernelAbstractions), with boolean shapes as flat CSG programs and BVHs for the
      daughters of volumes. `xray` images on CPUs and GPUs.
    - Single precision (`Float32`) support, needed by Apple GPUs.
    - Shape benchmarks with the methodology of VecGeom (`benchmark_shape`, `vecgeom_suite`) and
      comparison with VecGeom.
    - Documentation with Documenter and Literate examples.
- Changes:
    - `snake_case` API (`distance_to_in`, `read_gdml`, `world_volume`, ...) and own point and vector
      types (`Point3`, `Vector3`).
    - Drawing functions moved to a Makie package extension.
- Fixes:
    - Many fixes of the shape functions at surfaces and for grazing rays.
    - Boolean shapes: `inside` and `normal` at the points on the surface of both operands (seams of
      unions, edges where the cut of a subtraction meets the surface), and the quick bounding box
      checks for the points on the surface of an operand.
