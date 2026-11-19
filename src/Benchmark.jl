#---Benchmarks of the shape functions---------------------------------------------------------------
#   The methodology follows the shape benchmarks of VecGeom (test/VecGeomTest/Benchmarker.cpp and
#   VecGeom/volumes/utilities/VolumeUtilities.h), so that the timings can be compared:
#   - the shape is placed in a world box; points are sampled uniformly in the box ((1 - 2u) × half
#     length per axis) and directions uniformly in a cube and normalized (not isotropic);
#   - inside and contains: points in the world, a fraction `inside_bias` (0.5) of them in the shape;
#   - distance_to_in and safety_to_in: points in the world outside the shape, with directions such
#     that a fraction `toin_bias` (0.8) of the rays hit the shape;
#   - distance_to_out and safety_to_out: points inside the shape (in its own frame), random directions;
#   - each function is called for the `npoints` points, `nrep` times; the points in the world are
#     transformed to the frame of the shape in the timed loop, as for the placed volumes of VecGeom.
#   Two modes: `:shape` calls the functions of the shape directly (as the "Specialized" mode of
#   VecGeom), `:model` goes through the dispatch of a flat GeometryModel (as the navigator does).
#---------------------------------------------------------------------------------------------------
using Printf, Random
export ShapeBenchmark, benchmark_shape, vecgeom_suite

"""
    ShapeBenchmark{T}(name, shape, world, placement, npoints, nrep, vecgeom)

A benchmark case: a shape placed with `placement` (a `Transformation3D` as given to
`place_daughter!`: `Transformation3D(x, y, z)` puts the shape at (x, y, z)) in a world box of half
lengths `world`, timed for `npoints` points and `nrep` repetitions.
`vecgeom` is the VecGeom benchmark program (and options) with the same setup.
"""
struct ShapeBenchmark{T<:AbstractFloat}
    name::String
    shape::AbstractShape{T}
    world::Vector3{T}
    placement::Transformation3D{T}
    npoints::Int
    nrep::Int
    vecgeom::String
end

#---Sampling as in VecGeom----------------------------------------------------------------------------
@inline _sample_point(rng, c::StaticVector{3,T}, h::StaticVector{3,T}) where T =
    Point3{T}(c[1] + h[1] * (1 - 2rand(rng, T)), c[2] + h[2] * (1 - 2rand(rng, T)), c[3] + h[3] * (1 - 2rand(rng, T)))
@inline function _sample_direction(rng, ::Type{T}) where T
    d = Vector3{T}(1 - 2rand(rng, T), 1 - 2rand(rng, T), 1 - 2rand(rng, T))
    return d / sqrt(d ⋅ d)
end
_contained(shape, p) = inside(shape, p) != kOutside

# A point inside the shape, in its own frame (rejection in its bounding box, as FillRandomPoints)
function _point_inside(rng, shape::AbstractShape{T}) where T
    lo, hi = extent(shape)
    c, h = (lo + hi) / 2, (hi - lo) / 2
    while true
        p = _sample_point(rng, c, h)
        inside(shape, p) == kInside && return p
    end
end

# Points for inside/contains (world frame), as FillContainedPoints(world, bias, placed=true).
# Note: as in VecGeom, the points added to reach the bias are sampled in the frame of the shape and
# used as world points (the transformation applied is the one of the world, the identity).
function _inside_points(rng, shape, tr, world::Vector3{T}, n, bias) where T
    o = zero(Point3{T})
    pts = [_sample_point(rng, o, world) for _ in 1:n]
    isin = [_contained(shape, tr * p) for p in pts]
    count = sum(isin)
    i = 1
    while count / n > bias                       # replace contained points
        while !isin[i]; i += 1; end
        p = _sample_point(rng, o, world)
        while _contained(shape, tr * p)
            p = _sample_point(rng, o, world)
        end
        pts[i], isin[i] = p, false
        count -= 1
    end
    i = 1
    while count / n < bias                       # add contained points
        while isin[i]; i += 1; end
        pts[i], isin[i] = _point_inside(rng, shape), true
        count += 1
    end
    return pts
end

# Points for distance_to_in (world frame, not in the shape) and directions hitting the shape for a
# fraction `bias` of them, as FillUncontainedPoints and FillBiasedDirections
function _toin_points(rng, shape, tr, world::Vector3{T}, n, bias) where T
    o = zero(Point3{T})
    pts = Vector{Point3{T}}(undef, n)
    for i in 1:n
        p = _sample_point(rng, o, world)
        while _contained(shape, tr * p)
            p = _sample_point(rng, o, world)
        end
        pts[i] = p
    end
    hits(p, d) = distance_to_in(shape, tr * p, tr * d) < Inf
    dirs = [_sample_direction(rng, T) for _ in 1:n]
    hit = [hits(pts[i], dirs[i]) for i in 1:n]
    nhit = sum(hit)
    while nhit / n > bias                        # remove hits: invert or resample the direction
        k = rand(rng, 1:n)
        tries = 0
        while hit[k] && tries < 100
            tries += 1
            dirs[k] = isodd(tries) ? -dirs[k] : _sample_direction(rng, T)
            hit[k] = hits(pts[k], dirs[k])
            hit[k] || (nhit -= 1)
        end
    end
    itr = inv(tr)
    while nhit / n < bias                        # add hits: aim at a point of the shape
        k = rand(rng, 1:n)
        hit[k] && continue
        # VecGeom aims at a point on the surface (SamplePointOnSurface); here, a point inside
        d = itr * _point_inside(rng, shape) - pts[k]
        dirs[k] = d / sqrt(d ⋅ d)
        hit[k] = hits(pts[k], dirs[k])
        hit[k] && (nhit += 1)
    end
    return pts, dirs
end

# Points for distance_to_out (frame of the shape, inside it) and random directions
_toout_points(rng, shape::AbstractShape{T}, n) where T =
    [_point_inside(rng, shape) for _ in 1:n], [_sample_direction(rng, T) for _ in 1:n]

#---Timing---------------------------------------------------------------------------------------------
# The function `f` is a type parameter, so that the loops are specialized (no dynamic dispatch). Each
# call goes through a function that is not inlined: otherwise the compiler vectorizes the loop over
# the points for the simplest functions, which measures something else than a call (VecGeom makes a
# virtual call per point, and navigation calls the shapes one at a time).
@noinline _call(f::F, p) where F = f(p)
@noinline _call(f::F, p, d) where F = f(p, d)
function _time(f::F, out, nrep, pts) where F
    t0 = time_ns()
    for _ in 1:nrep, i in eachindex(pts)
        @inbounds out[i] = _call(f, pts[i])
    end
    return (time_ns() - t0) * 1e-9
end
function _time(f::F, out, nrep, pts, dirs) where F
    t0 = time_ns()
    for _ in 1:nrep, i in eachindex(pts)
        @inbounds out[i] = _call(f, pts[i], dirs[i])
    end
    return (time_ns() - t0) * 1e-9
end

# The shape functions in the two modes. Callable structs (not closures assigned in branches, which
# would not be type stable) so that the timed loops are fully specialized.
struct _ShapeCall{S}
    shape::S
end
@inline (c::_ShapeCall)(f::F, args...) where F = f(c.shape, args...)
struct _ModelCall{M}
    model::M
    ref::ShapeRef
end
@inline (c::_ModelCall)(f::F, args...) where F = shape_call(f, c.model, c.ref, args...)

_functions(call::C, tr) where C =
    (inside          = p -> call(inside, tr * p),
     contains        = p -> call(inside, tr * p) == kInside,
     distance_to_in  = (p, d) -> call(distance_to_in, tr * p, tr * d),
     safety_to_in    = p -> call(safety_to_in, tr * p),
     distance_to_out = (p, d) -> call(distance_to_out, p, d),
     safety_to_out   = p -> call(safety_to_out, p))

function _functions(shape, tr, mode::Symbol)
    mode == :shape && return _functions(_ShapeCall(shape), tr)
    model = close_geometry(_benchmark_world(shape, tr))
    return _functions(_ModelCall(model, model.volumes[model.placements[1].volume].shape), tr)
end

function _benchmark_world(shape::AbstractShape{T}, tr) where T
    lo, hi = extent(shape)
    h = 4 * max.(abs.(lo), abs.(hi)) .+ 1
    world = Volume{T}("world", Box{T}(h...), Material{T}("vacuum"; density=0))
    place_daughter!(world, tr, Volume{T}("shape", shape, Material{T}("material"; density=1)))
    return world
end

"""
    benchmark_shape(shape; world, placement, npoints=1024, nrep=4, inside_bias=0.5, toin_bias=0.8,
                    modes=(:shape, :model), seed=1, verbose=true)
    benchmark_shape(b::ShapeBenchmark; kwargs...)

Time the shape functions as the shape benchmarks of VecGeom: `inside`, `contains`, `distance_to_in`,
`safety_to_in`, `distance_to_out` and `safety_to_out` (see the header of Benchmark.jl). `world` are
the half lengths of the world box (default: 4 times the extent of the shape) and `placement` the
placement of the shape in it (as for `place_daughter!`). Returns a vector of named tuples (case, T, mode, op, npoints, nrep,
seconds, ns), with the total time and the time per call in ns.
"""
function benchmark_shape(shape::AbstractShape{T}; name::AbstractString=string(nameof(typeof(shape))),
                         world=nothing, placement::Transformation3D{T}=one(Transformation3D{T}),
                         npoints::Integer=1024, nrep::Integer=4, inside_bias::Real=0.5, toin_bias::Real=0.8,
                         modes=(:shape, :model), seed::Integer=1, verbose::Bool=true) where T
    if world === nothing
        lo, hi = extent(shape)
        world = 4 * max.(abs.(lo), abs.(hi))
    end
    world = Vector3{T}(world...)
    tr = placement                                # maps the world frame to the frame of the shape
    rng = MersenneTwister(seed)
    pin = _inside_points(rng, shape, tr, world, npoints, inside_bias)
    ptoin, dtoin = _toin_points(rng, shape, tr, world, npoints, toin_bias)
    ptoout, dtoout = _toout_points(rng, shape, npoints)
    rows = NamedTuple[]
    verbose && @printf "%s{%s}: %s, %d points × %d repetitions\n" name T string(shape) npoints nrep
    for mode in modes
        fs = _functions(shape, tr, mode)
        ints, bools, dists = Vector{Int}(undef, npoints), Vector{Bool}(undef, npoints), Vector{T}(undef, npoints)
        cases = ((:inside, fs.inside, ints, (pin,)), (:contains, fs.contains, bools, (pin,)),
                 (:distance_to_in, fs.distance_to_in, dists, (ptoin, dtoin)), (:safety_to_in, fs.safety_to_in, dists, (ptoin,)),
                 (:distance_to_out, fs.distance_to_out, dists, (ptoout, dtoout)), (:safety_to_out, fs.safety_to_out, dists, (ptoout,)))
        verbose && @printf "  %-6s" mode
        for (op, f, out, args) in cases
            _time(f, out, 1, args...)            # compile and warm up
            t = _time(f, out, nrep, args...)
            ns = 1e9 * t / (npoints * nrep)
            push!(rows, (case=name, T=string(T), mode=mode, op=op, npoints=npoints, nrep=nrep, seconds=t, ns=ns))
            verbose && @printf " %s %.1f ns" op ns
        end
        verbose && println()
    end
    return rows
end
benchmark_shape(b::ShapeBenchmark; kwargs...) =
    benchmark_shape(b.shape; name=b.name, world=b.world, placement=b.placement, npoints=b.npoints, nrep=b.nrep, kwargs...)

#---The benchmark programs of VecGeom (test/benchmark/*Benchmark.cpp) with their default setup--------
"""
    vecgeom_suite(T=Float64)

The shape benchmarks of VecGeom (test/benchmark/*Benchmark.cpp) for the shapes of Geom4hep, with the
same shape parameters, world box, placement, number of points and repetitions as the defaults of
the VecGeom programs. Not reproduced: the corner-defined Trapezoid (type 1) and the (r,z) Polycone
(type 4).
"""
function vecgeom_suite(::Type{T}=Float64) where T
    B(name, shape, world, placement, npoints, nrep, prog) = ShapeBenchmark{T}(name, shape, Vector3{T}(world...), placement, npoints, nrep, prog)
    tr(x, y, z) = Transformation3D{T}(x, y, z)
    deg = π / 180
    trapezoids = ((0, Trap{T}(10, 0, 0, 10, 10, 10, 0, 10, 10, 10, 0)),
                  (2, Trap{T}(15, 20deg, 5deg, 10, 7.5, 10, 10deg, 4, 2.5, 3.5, 10deg)),
                  (3, Trap{T}(15, 0, 0, 10, 7.5, 7.5, 0, 10, 3.5, 3.5, 0)),
                  (4, Trap{T}(15, 0, 0, 10, 7.5, 7.5, 0, 8, 3.5, 3.5, 0)),
                  (5, Trap{T}(15, 20deg, 5deg, 10, 7.5, 7.5, 10deg, 10, 7.5, 7.5, 10deg)))
    # VecGeom: Transformation3D(5, 2, 3, 15, 30, 45), Euler angles (ϕ, θ, ψ) in degrees
    traptr = Transformation3D{T}(5, 2, 3, Rotations.RotZXZ{T}(15deg, 30deg, 45deg))
    pcone(rmin, rmax, z) = Polycone{T}(rmin, rmax, z, 0, 2π)
    polycones = ((0, pcone([0.1, 0, 0, 0.2], [1, 2, 2, 1.5], [-1, -0.5, 0.5, 10])),
                 (1, pcone([0.2, 0], [1, 1.5], [-1, 3])),
                 (2, pcone([0, 0, 0, 0], [1, 2, 2, 1.5], [-1, -0.5, 0.5, 10])),
                 (3, pcone([0.1, 0.2, 0.2, 0.5], [1, 2, 2, 1.5], [-1, -0.5, 0.5, 10])))
    pgon(n, rmin, rmax, z, ϕ₀=0, Δϕ=2π) = Polyhedra{T}(n, rmin, rmax, z, ϕ₀, Δϕ)
    polyhedra = ((0, pgon(5, [0, 0, 0, 0, 0], [2, 3, 2, 3, 2], [-4, -2, 0, 2, 4])),
                 (1, pgon(5, [1, 0.75, 0.5, 0.75, 1], fill(1.5, 5), [-4, -1, 0, 1, 4])),
                 (2, pgon(5, [1, 0.75, 0.5, 0.75, 1], fill(1.5, 5), [-4, -1, 0, 1, 4], 15deg, 45deg)),
                 (3, pgon(5, [1, 0.75, 0.5, 0.75, 1], fill(1.5, 5), [-4, -1, 0, 1, 4], 15deg, 340deg)),
                 (4, pgon(6, fill(1, 17), [2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2], collect(-8:8))),
                 (5, pgon(5, [0, 1, 0.5, 1, 0], [1, 2, 2, 2.5, 1], [-2, -1, 1, 1, 2], 15deg, 340deg)))
    suite = [
        B("Box", Box{T}(1, 2, 3), (4, 8, 12), tr(0.1, 0, 0), 1024, 4, "BoxBenchmark"),
        B("Trd", Trd{T}(3, 4, 5, 6, 7), (16, 24, 28), tr(5, 5, 5), 50000, 1, "TrdBenchmark"),
        [B("Trap type $k", s, (20, 20, 20), traptr, 32768, 4, "TrapezoidBenchmark -type $k") for (k, s) in trapezoids]...,
        B("Tube", Tube{T}(0, 5, 10, 0, 2π), (20, 20, 40), tr(5, 5, 5), 10240, 1, "TubeBenchmark"),
        B("Cone", Cone{T}(5, 10, 7, 15, 10, 0, 2π), (100, 100, 100), tr(5, 0, 0), 10240, 1, "ConeBenchmark"),
        B("CutTube", CutTube{T}(3, 5, 10, 0, 2π/3, 3π/4, π/3, π/4, 2π/3), (20, 20, 40), tr(5, 5, 5), 10240, 1, "CutTubeBenchmark"),
        B("Sphere", Sphere{T}(15, 20, 0, 2π, 0, π), (80, 80, 80), tr(5, 5, 5), 1024, 10, "SphereBenchmark"),
        B("Orb", Sphere{T}(0, 3), (12, 12, 12), tr(5, 5, 5), 10240, 1, "OrbBenchmark"),
        [B("Polycone type $k", s, (10, 10, 35), tr(0, 0, 0), 1024, 4, "PolyconeBenchmark -type $k") for (k, s) in polycones]...,
        [B("Polyhedron type $k", s, (5, 5, 10), tr(0, 0, 0), 10000, 4, "PolyhedronBenchmark -type $k") for (k, s) in polyhedra]...,
        B("Torus", Torus{T}(1.2, 3.1, 5, 0, 2π), (16.2, 16.2, 16.2), tr(0, 0, 0), 1024, 4, "TorusBenchmark2"),
        B("MultiUnion", _random_boxes(T, 100), (20, 20, 20), tr(0, 0, 0), 1024, 4, "MultiUnionBenchmark"),
    ]
    return suite
end

# MultiUnion of `n` boxes with random sizes, positions and orientations, as MultiUnionBenchmark
function _random_boxes(::Type{T}, n; size=10, seed=1) where T
    rng = MersenneTwister(seed)
    sized = size * (0.5 / n)^(1 / 3)
    mu = MultiUnion{T}(Box{T}[], Transformation3D{T}[])
    u(a, b) = a + (b - a) * rand(rng)
    for _ in 1:n
        s = u(0.8sized, 1.2sized)
        pos = (u(-size, size), u(-size, size), u(-size, size))
        rot = Rotations.RotZXZ{T}(deg2rad(u(-180, 180)), deg2rad(u(-180, 180)), deg2rad(u(-180, 180)))
        add_component!(mu, Box{T}(s, s, s), Transformation3D{T}(pos..., rot))
    end
    return mu
end
