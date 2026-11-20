using Geom4hep, Geom4hep.Units
using Geom4hep.Rotations: RotZXZ
using LinearAlgebra
using CairoMakie: CairoMakie, Figure, Axis, LScene, heatmap!, Colorbar

const T = Float64;

g4transform(x, y, z, φ=0, θ=0, ψ=0) = Transformation3D{T}(x, y, z, RotZXZ{T}(φ, θ, ψ));

function roundcube(a, r)
    ca = a / 4
    cube  = Box{T}(ca, ca, ca)
    cyl   = Tube{T}(0, r, ca, 0, 2π)
    orb   = Sphere{T}(0, r)

    edge  = Box{T}(r, r, ca)
    vert  = Box{T}(r, r, r)
    frame = Box{T}(ca, ca, ca)

    # the material to remove along an edge (acyl) and at a corner (aorb)
    acyl1 = BooleanSubtraction(edge, cyl, g4transform(-r, -r, 0))
    acyl2 = BooleanSubtraction(edge, cyl, g4transform(-r, r, 0))
    acyl3 = BooleanSubtraction(edge, cyl, g4transform(r, r, 0))
    aorb  = BooleanSubtraction(vert, orb, g4transform(-r, -r, -r))

    # a corner of the cube with three rounded edges and a rounded vertex
    s1 = BooleanSubtraction(cube, acyl1, g4transform(ca, ca, 0))
    s2 = BooleanSubtraction(s1, acyl2, g4transform(ca, 0, ca, 0, π/2, 0))
    s3 = BooleanSubtraction(s2, acyl3, g4transform(0, ca, ca, 0, π/2, π/2))
    s4 = BooleanSubtraction(s3, aorb, g4transform(ca, ca, ca))

    # four corners around the frame make half of the cube (x > 0)
    u1 = BooleanUnion(frame, s4, g4transform(ca, ca, ca))
    u2 = BooleanUnion(u1, s4, g4transform(ca, ca, -ca, 0, π/2, 0))
    u3 = BooleanUnion(u2, s4, g4transform(ca, -ca, ca, 0, -π/2, 0))
    u4 = BooleanUnion(u3, s4, g4transform(ca, -ca, -ca, 0, π/2, π/2))

    # and the other half, rotated by π around z
    return BooleanUnion(u4, u4, g4transform(0, 0, 0, 0, 0, π))
end

a, r = 10cm, 1cm
rcube = roundcube(a, r);

extent(rcube)

function sdf(p, h, r)
    q = abs.(p) .- (h - r)
    return norm(max.(q, 0)) - r + min(maximum(q), 0)
end
points = [Point3{T}((rand(3) .- 0.5) .* 1.2a) for _ in 1:100_000]
mismatches = count(points) do p
    d = sdf(p, a/2, r)
    abs(d) > 1e-6 && (inside(rcube, p) == kOutside) != (d > 0)
end

b = a - 2r
exact = b^3 + 6b^2 * r + 3π * r^2 * b + 4π/3 * r^3
capacity(rcube) / exact

dir = normalize(Vector3{T}(-1, -1, -1))
p = Point3{T}(10cm, 10cm, 10cm)
distance_to_in(rcube, p, dir), norm(p) - (√3 * (a/2 - r) + r)

normal(rcube, Point3{T}(a/2 - r + r/√2, a/2 - r + r/√2, 1cm))

p_seam, p_face = Point3{T}(0, 0, 3cm), Point3{T}(a/2, 0, 0)
inside(rcube, p_seam) == kInside, inside(rcube, p_face) == kSurface, normal(rcube, p_face)

function surface_point(q, h, r)                     # projection on the surface, and the normal
    c = clamp.(q, -(h - r), h - r)
    d = q - c
    norm(d) > 0 && return c + r * normalize(d), normalize(d)
    i = argmax(abs.(q))
    n = Vector3{T}(ntuple(j -> j == i ? sign(q[i]) : 0, 3))
    return Point3{T}(ntuple(j -> j == i ? sign(q[i]) * h : q[j], 3)), n
end
seams = (0, a/4, -a/4)
bad = 0
for _ in 1:5000
    q = (rand(3) .- 0.5) .* 1.2a
    for _ in 1:rand(0:2)
        q[rand(1:3)] = rand(seams)                  # move to one or two of these planes
    end
    p, n = surface_point(Point3{T}(q...), a/2, r)
    global bad += inside(rcube, p) != kSurface || normal(rcube, p) ⋅ n < 0.9999
end
bad

draw_distance_to_in(rcube, 50_000)

xs = range(0, 0.6a, length=300)
fig = Figure(size=(900, 400))
for (i, z) in enumerate((0, a/2 - r/2))
    ax = Axis(fig[1, i], aspect=1, title="z = $(z) mm", xlabel="x (mm)", ylabel="y (mm)")
    heatmap!(ax, xs, xs, [Int(inside(rcube, Point3{T}(x, y, z))) for x in xs, y in xs], colormap=:grays)
end
fig

world = Volume{T}("world", Box{T}(40cm, 40cm, 15cm), Material{T}("air"; density=0.0012))
for (i, (radius, density)) in enumerate(((0.5cm, 2.7), (1.5cm, 7.87), (3cm, 11.35), (4cm, 1.0)))
    vol = Volume{T}("rcube$i", roundcube(20cm, radius), Material{T}("m$i"; density=density))
    x, y = (i - 1) % 2 == 0 ? -18cm : 18cm, (i - 1) ÷ 2 == 0 ? -18cm : 18cm
    place_daughter!(world, Transformation3D{T}(x, y, 0, 0, 0, π/8), vol)
end
model = close_geometry(world)
u, v, image, _ = xray(model, 3, 400)
fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1], aspect=1, xlabel="x (mm)", ylabel="y (mm)")
hm = heatmap!(ax, u, v, image, colormap=:grays)
Colorbar(fig[1, 2], hm, label="mass thickness (g/cm³ × mm)")
fig
