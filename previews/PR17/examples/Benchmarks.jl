using Geom4hep, Printf
using CairoMakie: CairoMakie, Figure, Axis, barplot!

rows = benchmark_shape(Torus{Float64}(1.2, 3.1, 5, 0, 2π); npoints=1024, nrep=10, verbose=false)
for r in rows
    @printf "%-6s %-16s %8.1f ns\n" r.mode r.op r.ns
end

suite = vecgeom_suite(Float64)
[(b.name, b.vecgeom) for b in suite][1:6]

results = Dict(b.name => benchmark_shape(b; verbose=false) for b in suite)
ops = (:inside, :distance_to_in, :distance_to_out, :safety_to_in, :safety_to_out)
@printf "%-18s" "case"; foreach(op -> @printf("%16s", op), ops); println()
for b in suite
    @printf "%-18s" b.name
    for op in ops
        @printf "%16.1f" only(r.ns for r in results[b.name] if r.mode == :shape && r.op == op)
    end
    println()
end

names = [b.name for b in suite]
times = [only(r.ns for r in results[n] if r.mode == :shape && r.op == :distance_to_in) for n in names]
fig = Figure(size=(900, 400))
ax = Axis(fig[1, 1], ylabel="ns per call", yscale=log10, xticks=(1:length(names), names),
          xticklabelrotation=π/3, title="distance_to_in")
barplot!(ax, 1:length(names), times)
fig
