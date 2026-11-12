#---Navigation in a flat GeometryModel-----------------------------------------------------------------
#   The state is an immutable isbits struct and the functions return new states, so that the same
#   code runs on the CPU and in GPU kernels. Daughters are visited with an explicit iterator (linear
#   scan or stackless traversal of the BVH), without recursion. The main functions are not inlined:
#   with a constant starting depth the compilers would fully unroll the loops over the levels of the
#   hierarchy, which makes GPU kernels explode in size.
#---------------------------------------------------------------------------------------------------

"""
    NavState{T}

Position of a point in the volume hierarchy of a `GeometryModel`: the current volume, the path of
placements from the world (`depth` levels) and the transformation from the global frame to the frame
of the current volume. `volume == 0` when the point is outside the world.
"""
struct NavState{T<:AbstractFloat}
    depth::Int32
    volume::Int32
    path::NTuple{MAX_NAV_DEPTH,Int32}
    tolocal::Transformation3D{T}
end

@inline _world_state(::Type{T}) where T = NavState{T}(0, 1, ntuple(_ -> Int32(0), Val(MAX_NAV_DEPTH)), one(Transformation3D{T}))
@inline _outside_state(::Type{T}) where T = NavState{T}(-1, 0, ntuple(_ -> Int32(0), Val(MAX_NAV_DEPTH)), one(Transformation3D{T}))

"`true` if the state corresponds to a point inside the world volume."
@inline in_world(s::NavState) = s.volume != 0

"Index of the material of the current volume."
@inline material_index(model::GeometryModel, s::NavState) = @inbounds model.volumes[s.volume].material

#---Moving along the hierarchy-------------------------------------------------------------------------
@inline function _push(model::GeometryModel{T}, s::NavState{T}, pidx::Int32) where T
    pl = @inbounds model.placements[pidx]
    NavState{T}(s.depth + Int32(1), pl.volume, Base.setindex(s.path, pidx, s.depth + 1), pl.transformation * s.tolocal)
end

@noinline function _pop(model::GeometryModel{T}, s::NavState{T}) where T
    s.depth <= 0 && return _outside_state(T)
    depth = s.depth - Int32(1)
    volume = Int32(1)
    tolocal = one(Transformation3D{T})
    for level in 1:depth                 # recompute the transformation from the world
        pl = @inbounds model.placements[s.path[level]]
        tolocal = pl.transformation * tolocal
        volume = pl.volume
    end
    NavState{T}(depth, volume, Base.setindex(s.path, Int32(0), s.depth), tolocal)
end

#---Iteration over the daughters of a volume (filtered by a test on their bounding boxes)-------------
struct _DaughterIterator
    node::Int32      # BVH: next node to visit (0 = done); linear: next daughter (1-based)
    j::Int32         # BVH: next entry of the current leaf
    jend::Int32      # BVH: last entry of the current leaf
end
@inline _daughters_start(vol::FlatVolume) = _DaughterIterator(vol.bvh == 0 ? Int32(1) : vol.bvh, Int32(1), Int32(0))

"Next daughter (placement index) whose bounding box passes `test`, or 0 when done."
function _next_daughter(model::GeometryModel, vol::FlatVolume, it::_DaughterIterator, test::F) where F
    if vol.bvh == 0
        k = it.node
        while k <= vol.ndaughters
            pidx = vol.first_daughter + k
            k += Int32(1)
            test(@inbounds(model.placements[pidx]).aabb) && return pidx, _DaughterIterator(k, it.j, it.jend)
        end
        return Int32(0), _DaughterIterator(k, it.j, it.jend)
    else
        node, j, jend = it.node, it.j, it.jend
        while true
            while j <= jend                                   # entries of the current leaf
                pidx = @inbounds model.bvhindices[j]
                j += Int32(1)
                test(@inbounds(model.placements[pidx]).aabb) && return pidx, _DaughterIterator(node, j, jend)
            end
            node == 0 && return Int32(0), _DaughterIterator(node, j, jend)
            n = @inbounds model.bvhnodes[node]
            if test(n.aabb)
                if n.count > 0                                # leaf: visit its entries, then escape
                    j, jend = n.first, n.first + n.count - Int32(1)
                    node = n.skip
                else                                          # internal: descend (left child follows)
                    node += Int32(1)
                end
            else
                node = n.skip
            end
        end
    end
end

#---Locating points------------------------------------------------------------------------------------
@inline _inside_volume(model::GeometryModel, vidx, lpoint) = shape_call(inside, model, @inbounds(model.volumes[vidx]).shape, lpoint)

# Go down the hierarchy to the deepest volume containing the point. With a direction `gdir` (the
# point moving along it), a daughter whose surface the point is on is also entered when the direction
# goes into it: a point pushed across a boundary shared by a daughter of the entered volume can stay
# within the surface tolerance of that daughter (e.g. crossing a cylinder at a grazing angle).
@noinline function _descend(model::GeometryModel{T}, s::NavState{T}, gpoint::Point3{T}, gdir=nothing) where T
    while s.depth < MAX_NAV_DEPTH
        vol = @inbounds model.volumes[s.volume]
        lpoint = s.tolocal * gpoint
        found = Int32(0)
        surface = Int32(0)                  # first daughter entered at its surface
        it = _daughters_start(vol)
        while true
            pidx, it = _next_daughter(model, vol, it, bb -> inside(bb, lpoint))
            pidx == 0 && break
            pl = @inbounds model.placements[pidx]
            dpoint = pl.transformation * lpoint
            pos = _inside_volume(model, pl.volume, dpoint)
            if pos == kInside
                found = pidx
                break
            elseif gdir !== nothing && pos == kSurface && surface == 0
                # entering: on the surface, distance_to_in is 0 or (within the tolerance) negative
                ddir = pl.transformation * (s.tolocal * gdir)
                d = shape_call(distance_to_in, model, @inbounds(model.volumes[pl.volume]).shape, dpoint, ddir)
                d <= 0 && (surface = pidx)
            end
        end
        found == 0 && (found = surface)
        found == 0 && return s
        s = _push(model, s, found)
    end
    return s
end

# Go up the hierarchy until a volume contains the point, then down again
@noinline function _relocate(model::GeometryModel{T}, s::NavState{T}, gpoint::Point3{T}, gdir=nothing) where T
    while in_world(s) && _inside_volume(model, s.volume, s.tolocal * gpoint) != kInside
        s = _pop(model, s)
    end
    in_world(s) || return s
    return _descend(model, s, gpoint, gdir)
end

"""
    locate(model::GeometryModel, point) -> NavState

Find the deepest volume containing the (global) `point`.
"""
@noinline function locate(model::GeometryModel{T}, gpoint::Point3{T}) where T
    _inside_volume(model, Int32(1), gpoint) == kOutside && return _outside_state(T)
    return _descend(model, _world_state(T), gpoint)
end

#---Computing steps-----------------------------------------------------------------------------------
# Smallest step that moves a point: the push tolerance, but at least a few ulps of its coordinates
# (relevant for Float32 far from the origin, where steps of push_tolerance can be lost in rounding)
@inline _min_step(p::Point3{T}) where T = max(push_tolerance(T), 8 * eps(max(abs(p[1]), abs(p[2]), abs(p[3]))))

"""
    compute_step(model::GeometryModel, state, point, dir, step_limit) -> (step, new_state)

Distance from the global `point` along `dir` to the next boundary (or `step_limit`), and the state
after crossing it. The new state is computed for the point pushed slightly beyond the boundary.
"""
@noinline function compute_step(model::GeometryModel{T}, s::NavState{T}, gpoint::Point3{T}, gdir::Vector3{T}, step_limit::T) where T
    in_world(s) || return T(Inf), s
    vol = @inbounds model.volumes[s.volume]
    lpoint, ldir = s.tolocal * gpoint, s.tolocal * gdir
    dout = shape_call(distance_to_out, model, vol.shape, lpoint, ldir)
    if dout < 0
        # The point is (numerically) outside the current volume: move by a small step and relocate,
        # so that the caller always makes progress
        step = _min_step(gpoint)
        return step, _relocate(model, s, gpoint + step * gdir, gdir)
    end
    step = min(dout, step_limit)
    hit = Int32(0)
    invdir = inv.(ldir)
    it = _daughters_start(vol)
    while true
        pidx, it = _next_daughter(model, vol, it, bb -> intersect(bb, lpoint, ldir, invdir))
        pidx == 0 && break
        pl = @inbounds model.placements[pidx]
        d = shape_call(distance_to_in, model, @inbounds(model.volumes[pl.volume]).shape, pl.transformation * lpoint, pl.transformation * ldir)
        # d == 0: on the surface of the daughter, typically the one just left (an adjacent volume is
        # entered by the relocation of the pushed point instead); d < 0: already inside, enter it
        d == 0 && continue
        d < 0 && (d = zero(T))
        if d < step
            step = d
            hit = pidx
        end
    end
    if hit == 0 && step < dout                     # limited by step_limit: still in the same volume
        return step, s
    end
    # Crossing a boundary: make sure the point moves forward (very small steps would leave the
    # caller on the boundary while the new state is computed for the pushed point)
    minstep = _min_step(gpoint)
    step = max(step, minstep)
    newpoint = gpoint + (step + minstep) * gdir
    if hit != 0                                    # entering a daughter
        return step, _descend(model, _push(model, s, hit), newpoint, gdir)
    else                                           # leaving the current volume
        return step, _relocate(model, _pop(model, s), newpoint, gdir)
    end
end
