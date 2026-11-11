#---X-ray image of a geometry (benchmark of the navigation, as the XRay benchmark of VecGeom)-----------
#   Rays parallel to one of the axes cross the world volume; each pixel accumulates the mass thickness
#   (density × length) of the traversed volumes. The same per-pixel function runs in a
#   KernelAbstractions kernel on GPUs, or in a loop over Julia threads on the CPU.
#
#   Note: the CPU backend of KernelAbstractions ≥ 0.10 (POCL) inlines all functions, which makes the
#   kernels explode in size for geometries with nested boolean shapes; Julia threads are used instead.

@inline function _xray_pixel(model::GeometryModel{T}, origin, du, dv, dir, maxsteps, i, j) where T
    point = origin + T(i - 1) * du + T(j - 1) * dv
    state = locate(model, point)
    mass = zero(T)
    n = Int32(0)
    while in_world(state) && n < maxsteps
        density = @inbounds model.materials[material_index(model, state)].density
        step, state = compute_step(model, state, point, dir, T(Inf))
        point = point + step * dir
        mass += step * density
        n += Int32(1)
    end
    return mass, n
end

@kernel function _xray_kernel!(image, nsteps, model, origin, du, dv, dir, maxsteps, j0)
    i, jj = @index(Global, NTuple)
    j = jj + j0                       # columns are processed in batches (see `batch` in `xray`)
    mass, n = _xray_pixel(model, origin, du, dv, dir, maxsteps, i, j)
    @inbounds image[i, j] = mass
    @inbounds nsteps[i, j] = n
end

"""
    xray(model::GeometryModel, axis, npixels; backend=CPU(), maxsteps=10_000, launch_time=0.5)

X-ray image of the world volume of `model` with rays parallel to `axis` (1, 2 or 3 for x, y, z):
the mass thickness (density × length) integrated along each ray. `npixels` is the number of pixels
along the first transverse axis. Returns `(u, v, image, nsteps)` with the pixel centres along the two
transverse axes, the image and the number of navigation steps per pixel.

`backend` is a KernelAbstractions backend: `CPU()` (Julia threads), or a GPU backend such as
`MetalBackend()` or `CUDABackend()` (the model is copied to the device). On GPUs the image is
computed in several kernel launches, each one taking about `launch_time` seconds: a long-running
kernel on a GPU that also drives a display is aborted by the operating system (e.g. macOS error
"Impacting Interactivity").
"""
function xray(model::GeometryModel{T}, axis::Integer, npixels::Integer; backend=CPU(), maxsteps::Integer=10_000,
              launch_time::Real=0.5) where T
    lower, upper = shape_call(extent, model, model.volumes[1].shape)
    iu, iv = axis == 1 ? (2, 3) : axis == 2 ? (1, 3) : (1, 2)
    nu = npixels
    nv = max(1, round(Int, nu * (upper[iv] - lower[iv]) / (upper[iu] - lower[iu])))
    δu, δv = (upper[iu] - lower[iu]) / nu, (upper[iv] - lower[iv]) / nv
    unit(k) = Vector3{T}(k == 1, k == 2, k == 3)
    origin = lower + push_tolerance(T) * unit(axis) + (δu / 2) * unit(iu) + (δv / 2) * unit(iv)
    du, dv, dir = δu * unit(iu), δv * unit(iv), unit(axis)
    if backend isa CPU
        image = zeros(T, nu, nv)
        nsteps = zeros(Int32, nu, nv)
        Threads.@threads for j in 1:nv
            for i in 1:nu
                image[i, j], nsteps[i, j] = _xray_pixel(model, origin, du, dv, dir, Int32(maxsteps), i, j)
            end
        end
    else
        dimage = KernelAbstractions.zeros(backend, T, nu, nv)
        dnsteps = KernelAbstractions.zeros(backend, Int32, nu, nv)
        dmodel = adapt(backend, model)
        kernel! = _xray_kernel!(backend)
        # Columns are processed in batches; the size of the batches adapts to take ~launch_time each
        ncols, j0 = max(1, 4096 ÷ nu), 0
        while j0 < nv
            nj = min(ncols, nv - j0)
            t = @elapsed begin
                kernel!(dimage, dnsteps, dmodel, origin, du, dv, dir, Int32(maxsteps), j0; ndrange=(nu, nj))
                KernelAbstractions.synchronize(backend)
            end
            j0 += nj
            ncols = clamp(round(Int, nj * launch_time / max(t, 1e-3)), 1, 2nj)   # grow at most ×2
        end
        image, nsteps = Array(dimage), Array(dnsteps)
    end
    u = range(lower[iu] + δu / 2, step=δu, length=nu)
    v = range(lower[iv] + δv / 2, step=δv, length=nv)
    return u, v, image, nsteps
end
