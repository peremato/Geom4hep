#  Benchmarks of the shape functions of Geom4hep, with the methodology and setups of the shape
#  benchmarks of VecGeom (see `vecgeom_suite` and src/Benchmark.jl), optionally compared with VecGeom.
#
#  Usage (from the repository root):
#     julia --project=examples examples/Benchmarks.jl [options] [name filter]
#  Options:
#     --vecgeom DIR   run also the VecGeom benchmark programs found in DIR (a VecGeom build configured
#                     with -DBUILD_TESTING=ON -DVECGEOM_TEST_BENCHMARK=ON: <build>/test/benchmark),
#                     with the same numbers of points and repetitions
#     --calls N       minimum number of calls of each function (default 1000000): the number of
#                     repetitions of each case is increased to reach it (the timings of VecGeom's
#                     defaults are too short to be precise)
#     --float32       Geom4hep in single precision (VecGeom in the precision of its build)
#     --csv FILE      file for all the results (default: benchmarks.csv)
#  Example:
#     julia --project=examples examples/Benchmarks.jl --vecgeom ~/vecgeom-build/test/benchmark Polyhedron

using Geom4hep, Printf

#---Arguments------------------------------------------------------------------------------------------
function option(name, default)
    i = findfirst(==(name), ARGS)
    i === nothing && return default
    v = ARGS[i+1]
    deleteat!(ARGS, i:i+1)
    return v
end
vgdir = option("--vecgeom", nothing)
calls = parse(Int, option("--calls", "1000000"))
csvfile = option("--csv", "benchmarks.csv")
T = "--float32" in ARGS ? Float32 : Float64
filter!(!=("--float32"), ARGS)
pattern = isempty(ARGS) ? "" : ARGS[1]

#---VecGeom--------------------------------------------------------------------------------------------
const vgops = (inside="Inside", contains="Contains", distance_to_in="DistanceToIn", safety_to_in="SafetyToIn",
               distance_to_out="DistanceToOut", safety_to_out="SafetyToOut")

# Run a VecGeom benchmark program and return the total times (s) of its Specialized and Unspecialized modes
function run_vecgeom(dir, program, npoints, nrep)
    prog, opts... = split(program)
    exe = joinpath(dir, prog)
    isfile(exe) || (@warn "VecGeom program $exe not found"; return nothing)
    out = mktempdir() do tmp                     # some programs write files in the working directory
        read(pipeline(Cmd(`$exe $opts -npoints $npoints -nrep $nrep`; dir=tmp); stderr=devnull), String)
    end
    times = Dict{Tuple{String,String},Float64}()
    for line in split(out, '\n')
        m = match(r"^(Specialized|Unspecialized)\s+-", line)
        m === nothing && continue
        for (op, t) in eachmatch(r"(\w+): ([0-9.eE+-]+)s", line)
            times[(m[1], op)] = parse(Float64, t)
        end
    end
    return times
end

#---Run------------------------------------------------------------------------------------------------
suite = filter(b -> occursin(pattern, b.name), vecgeom_suite(T))
@printf "Geom4hep shape benchmarks (%s)%s, at least %d calls per function; times in ns per call\n" T (vgdir === nothing ? "" : ", compared with VecGeom") calls
@printf "%-20s %-16s %10s %10s" "case" "function" "Geom4hep" "(model)"
vgdir === nothing || @printf " %10s %10s %8s" "VecGeom" "(unspec.)" "ratio"
println()
results = []
for b in suite
    nrep = max(b.nrep, cld(calls, b.npoints))
    rows = benchmark_shape(b; nrep, verbose=false)
    vg = vgdir === nothing ? nothing : run_vecgeom(vgdir, b.vecgeom, b.npoints, nrep)
    for (k, op) in enumerate(keys(vgops))
        shape = only(r for r in rows if r.mode == :shape && r.op == op)
        model = only(r for r in rows if r.mode == :model && r.op == op)
        @printf "%-20s %-16s %10.1f %10.1f" (k == 1 ? b.name : "") op shape.ns model.ns
        row = (case=b.name, T=string(T), npoints=b.npoints, nrep=nrep, op=string(op), geom4hep_ns=shape.ns, geom4hep_model_ns=model.ns,
               vecgeom_ns=NaN, vecgeom_unspecialized_ns=NaN)
        if vg !== nothing
            n = b.npoints * nrep
            s = 1e9 * get(vg, ("Specialized", vgops[op]), NaN) / n
            u = 1e9 * get(vg, ("Unspecialized", vgops[op]), NaN) / n
            @printf " %10.1f %10.1f %8.2f" s u shape.ns / s
            row = merge(row, (vecgeom_ns=s, vecgeom_unspecialized_ns=u))
        end
        println()
        push!(results, row)
    end
end
vgdir === nothing || println("ratio = Geom4hep / VecGeom (Specialized): below 1 Geom4hep is faster")

open(csvfile, "w") do io
    println(io, join(keys(first(results)), ","))
    for r in results
        println(io, join((v isa AbstractFloat ? @sprintf("%.3f", v) : string(v) for v in values(r)), ","))
    end
end
println("results written to $csvfile")
