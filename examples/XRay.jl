#  X-ray images of a detector: the mass thickness (density × length) integrated along rays parallel to
#  the axes, computed with the navigator of a flat GeometryModel, on the CPU or on a GPU.
#
#  Usage (from the repository root):
#     julia --project=examples -t auto examples/XRay.jl [cpu|metal|cuda] [Float32|Float64]
#  The GPU packages are not part of the examples environment: add Metal (Apple) or CUDA (NVIDIA)
#  to it first, e.g. `julia --project=examples -e 'using Pkg; Pkg.add("Metal")'`.

using Geom4hep, Printf, GLMakie

device = length(ARGS) > 0 ? ARGS[1] : "cpu"
T      = length(ARGS) > 1 ? eval(Symbol(ARGS[2])) : (device == "cpu" ? Float64 : Float32)

backend = if device == "metal"
    @eval using Metal
    MetalBackend()                         # Metal GPUs only support Float32
elseif device == "cuda"
    @eval using CUDA
    CUDABackend()
else
    Geom4hep.CPU()
end

#-----build the detector-----------------------------------------------------------------------------
full = read_gdml(joinpath(@__DIR__, "cms2018.gdml"), T)
world = world_volume(full[1,7,1])          # a muon barrel sector; use `full` for the complete detector
model = close_geometry(world)
println(model)

#-----generate the images----------------------------------------------------------------------------
xray(model, 1, 8; backend)                 # compile
images = map(1:3) do axis
    t = @elapsed u, v, img, nsteps = xray(model, axis, 1000; backend)
    @printf "axis %d: %d×%d pixels, %.3g steps in %.2f s on %s\n" axis size(img)... sum(nsteps) t device
    (u, v, img)
end

#-----plot the results-------------------------------------------------------------------------------
limits = (0, maximum(maximum(img) for (_, _, img) in images))
fig = Figure(size = (1000, 1000))
for (axis, (u, v, img)) in enumerate(images)
    pos = ((1, 1), (1, 2), (2, 1))[axis]
    heatmap!(Axis(fig[pos...], title = "rays along $("xyz"[axis])", aspect = DataAspect()), u, v, img,
             colormap = :grayC, colorrange = limits)
end
save("cms2018-xray.png", fig)
println("image saved in cms2018-xray.png")
