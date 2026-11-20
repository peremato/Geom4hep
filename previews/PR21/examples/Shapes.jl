using Geom4hep, Geom4hep.Units
using LinearAlgebra
using CairoMakie: CairoMakie, Figure, LScene

const T = Float64;

box    = Box{T}(10cm, 15cm, 20cm)
tube   = Tube{T}(5cm, 10cm, 20cm, 0, 270deg)
cone   = Cone{T}(2cm, 5cm, 5cm, 10cm, 20cm, 0, 2π)
trd    = Trd{T}(10cm, 5cm, 15cm, 8cm, 20cm)
trap   = Trap{T}(20cm, 10deg, 30deg, 10cm, 8cm, 12cm, 5deg, 6cm, 5cm, 7cm, 5deg)
sphere = Sphere{T}(8cm, 10cm, 0, 270deg, 30deg, 120deg)
torus  = Torus{T}(2cm, 4cm, 12cm, 0, 270deg)
pcone  = Polycone{T}([0, 5cm, 5cm, 0], [8cm, 10cm, 10cm, 6cm], [-20cm, -10cm, 10cm, 20cm], 0, 2π)
phedra = Polyhedra{T}(6, [3cm, 3cm, 3cm], [6cm, 10cm, 6cm], [-20cm, 0, 20cm], 0, 2π)
ctube  = CutTube{T}(5cm, 10cm, 20cm, 0, 2π, Vector3{T}(0, -0.5, -1), Vector3{T}(0.5, 0, 1));

torus

shapes = [box, tube, cone, trd, trap, sphere, torus, pcone, phedra, ctube]
fig = Figure(size=(1000, 450))
for (i, s) in enumerate(shapes)
    scene = LScene(fig[(i-1) ÷ 5 + 1, (i-1) % 5 + 1], show_axis=false)
    draw!(scene, s)
end
fig

p_in, p_out = Point3{T}(0, 0, 0), Point3{T}(30cm, 0, 0)
inside(box, p_in) == kInside, inside(box, p_out) == kOutside, inside(box, Point3{T}(10cm, 0, 0)) == kSurface

dir = Vector3{T}(-1, 0, 0)
distance_to_in(box, p_out, dir), distance_to_out(box, p_in, dir)

distance_to_in(box, p_out, Vector3{T}(0, 1, 0))

safety_to_in(box, p_out), safety_to_out(box, Point3{T}(5cm, 0, 0))

normal(sphere, Point3{T}(0, 10cm, 0))

extent(torus)

capacity(box) ≈ 8 * 10cm * 15cm * 20cm, surface(box)

draw_distance_to_out(phedra, 50_000)

hole  = Tube{T}(0, 4cm, 30cm, 0, 2π)
frame = BooleanSubtraction(box, hole)
frame = BooleanSubtraction(frame, hole, Transformation3D{T}(0, 0, 0, π/2, 0, 0))
frame = BooleanUnion(frame, Sphere{T}(0, 6cm), Transformation3D{T}(0, 0, 20cm))
inside(frame, Point3{T}(0, 0, 0)) == kOutside, inside(frame, Point3{T}(0, 0, 24cm)) == kInside

draw_distance_to_in(frame, 50_000)

box32 = Box{Float32}(10cm, 15cm, 20cm)
distance_to_in(box32, Point3{Float32}(30cm, 0, 0), Vector3{Float32}(-1, 0, 0)), tolerance(Float32)
