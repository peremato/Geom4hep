#---Build a flat GeometryModel from a hierarchy of volumes----------------------------------------------
const MAX_NAV_DEPTH = 24     # maximum depth of the volume hierarchy supported by the navigator

mutable struct _ModelBuilder{T}
    kinds::Dict{DataType,Int}             # concrete flat shape type => kind
    arrays::Vector{Vector}                # one array per kind
    shaperefs::IdDict{Any,ShapeRef}       # tree shape => reference (shared shapes are stored once)
    depth::Dict{ShapeRef,Int}             # boolean nesting level of each shape
    volumes::Vector{FlatVolume}
    volidx::IdDict{Any,Int32}
    placements::Vector{FlatPlacement{T}}
    materials::Vector{FlatMaterial{T}}
    matidx::IdDict{Any,Int32}
    bvhnodes::Vector{BVHNode{T}}
    bvhindices::Vector{Int32}
    cones::Vector{Cone{T}}
    zoffsets::Vector{T}
    csg_components::Vector{CSGComponent{T}}
    csg_program::Vector{Int32}
    bvh_threshold::Int
end

function _ModelBuilder{T}(bvh_threshold) where T
    _ModelBuilder{T}(Dict{DataType,Int}(), Vector[], IdDict{Any,ShapeRef}(), Dict{ShapeRef,Int}(),
                     FlatVolume[], IdDict{Any,Int32}(), FlatPlacement{T}[], FlatMaterial{T}[], IdDict{Any,Int32}(),
                     BVHNode{T}[], Int32[], Cone{T}[], T[], CSGComponent{T}[], Int32[], bvh_threshold)
end

function _store!(b::_ModelBuilder, flat, depth::Int)
    S = typeof(flat)
    isbitstype(S) || error("Shape type $S is not isbits and cannot be stored in a GeometryModel")
    kind = get!(b.kinds, S) do
        push!(b.arrays, S[])
        length(b.arrays)
    end
    push!(b.arrays[kind], flat)
    ref = ShapeRef(kind, length(b.arrays[kind]))
    b.depth[ref] = depth
    return ref
end

_aabb(shape, t::Transformation3D) = AABB(transform_extent(extent(shape), t)...)
function _tol_aabb(shape, t::Transformation3D{T}) where T
    tol = Vector3{T}(tolerance(T), tolerance(T), tolerance(T))
    lo, hi = transform_extent(extent(shape), t)
    AABB{T}(lo - tol, hi + tol)
end

#---Conversion of shapes---------------------------------------------------------------------------------
function _add_shape!(b::_ModelBuilder{T}, shape) where T
    haskey(b.shaperefs, shape) && return b.shaperefs[shape]
    ref = _convert_shape!(b, shape)
    b.shaperefs[shape] = ref
    return ref
end

_convert_shape!(b::_ModelBuilder{T}, shape::Union{Box{T},Trd{T},Trap{T},Tube{T},Cone{T},CutTube{T}}) where T = _store!(b, shape, 0)
_convert_shape!(b::_ModelBuilder{T}, shape::NoShape{T}) where T = _store!(b, FlatNoShape{T}(), 0)

function _convert_shape!(b::_ModelBuilder{T}, shape::Polycone{T}) where T
    first = length(b.cones)
    append!(b.cones, shape.sections)
    append!(b.zoffsets, shape.zᵢ)
    _store!(b, FlatPolycone{T}(first, length(shape.sections), shape.ϕ₀, shape.Δϕ), 0)
end

# Boolean shapes become CSG programs over primitive components
const _BooleanLike{T} = Union{BooleanShape{T}, MultiUnion{T}}

# Emit the program of `shape`, placed with `t` (CSG frame -> shape frame). Returns the stack depth used.
function _emit_csg!(b::_ModelBuilder{T}, comps, prog, shape, t::Transformation3D{T}) where T
    if shape isa BooleanShape
        op = shape isa BooleanUnion ? CSG_UNION : shape isa BooleanSubtraction ? CSG_SUBTRACTION : CSG_INTERSECTION
        dl = _emit_csg!(b, comps, prog, shape.left, t)
        dr = _emit_csg!(b, comps, prog, shape.right, shape.transformation * t)
        push!(prog, op)
        return max(dl, dr + 1)
    elseif shape isa MultiUnion
        depth = 0
        for i in 1:ncomponents(shape)
            d = _emit_csg!(b, comps, prog, shape.shapes[i], shape.transformations[i] * t)
            depth = max(depth, d + (i > 1))
            i > 1 && push!(prog, CSG_UNION)
        end
        return depth
    else
        push!(comps, CSGComponent{T}(_add_shape!(b, shape), t, _tol_aabb(shape, t)))
        push!(prog, Int32(length(comps)))
        return 1
    end
end

function _convert_shape!(b::_ModelBuilder{T}, shape::_BooleanLike{T}) where T
    comps, prog = CSGComponent{T}[], Int32[]
    depth = _emit_csg!(b, comps, prog, shape, one(Transformation3D{T}))
    depth <= CSG_MAX_STACK || error("Boolean shape too complex: the CSG program needs a stack of $depth > $CSG_MAX_STACK")
    lo, hi = extent(shape)
    csg = FlatCSG{T}(length(b.csg_components), length(comps), length(b.csg_program), length(prog), AABB{T}(lo, hi))
    append!(b.csg_components, comps)
    append!(b.csg_program, prog)
    _store!(b, csg, 1)
end

_convert_shape!(b::_ModelBuilder, shape) = error("Shape $(typeof(shape)) is not supported in a GeometryModel")

#---Materials, volumes and BVHs--------------------------------------------------------------------------
function _add_material!(b::_ModelBuilder{T}, mat) where T
    get!(b.matidx, mat) do
        push!(b.materials, FlatMaterial{T}(mat.density, mat.temperature, mat.Amass))
        Int32(length(b.materials))
    end
end

function _add_volume!(b::_ModelBuilder{T}, vol::Volume{T}, depth::Int) where T
    depth < MAX_NAV_DEPTH || error("The volume hierarchy is deeper than $MAX_NAV_DEPTH levels")
    haskey(b.volidx, vol) && return b.volidx[vol]
    vol.shape isa Aggregate && error("Volume $(vol.label) is an assembly: assemblies can only be placed")
    push!(b.volumes, FlatVolume(ShapeRef(0, 0), 0, 0, 0, 0))     # reserve the slot (preorder numbering)
    idx = Int32(length(b.volumes))
    b.volidx[vol] = idx
    shape = _add_shape!(b, vol.shape)
    material = _add_material!(b, vol.material)
    # daughters are stored contiguously
    first = Int32(length(b.placements))
    n = length(vol.daughters)
    for d in vol.daughters
        push!(b.placements, FlatPlacement{T}(0, d.transformation, d.aabb))
    end
    bvh = n > b.bvh_threshold ? _add_bvh!(b, vol, first) : Int32(0)
    for (i, d) in enumerate(vol.daughters)
        dv = _add_volume!(b, d.volume, depth + 1)
        p = b.placements[first + i]
        b.placements[first + i] = FlatPlacement{T}(dv, p.transformation, p.aabb)
    end
    b.volumes[idx] = FlatVolume(shape, material, first, n, bvh)
    return idx
end

# Flatten the (recursive) BVH of the daughters in depth-first order with escape indices: an internal
# node is followed by its left subtree and then its right subtree; `skip` is the node to continue
# with when the subtree is not visited (or after a leaf).
function _add_bvh!(b::_ModelBuilder{T}, vol::Volume{T}, first::Int32) where T
    root = Int32(length(b.bvhnodes) + 1)
    PENDING = Int32(-1)        # escape index not yet known
    function addnode!(node::BVH, skip::Int32)
        if isempty(node.children)
            push!(b.bvhnodes, BVHNode{T}(node.aabb, Int32(length(b.bvhindices) + 1), Int32(length(node.indices)), skip))
            append!(b.bvhindices, (first + Int32(i) for i in node.indices))
        else
            push!(b.bvhnodes, BVHNode{T}(node.aabb, 0, 0, skip))
            lstart = length(b.bvhnodes) + 1
            addnode!(node.children[1], PENDING)
            rstart = Int32(length(b.bvhnodes) + 1)
            for k in lstart:length(b.bvhnodes)            # the left subtree escapes to the right one
                n = b.bvhnodes[k]
                n.skip == PENDING && (b.bvhnodes[k] = BVHNode{T}(n.aabb, n.first, n.count, rstart))
            end
            addnode!(node.children[2], skip)
        end
    end
    addnode!(build_bvh(vol.daughters), Int32(0))
    return root
end

"""
    close_geometry(world::Volume; bvh_threshold=4)

Build a flat [`GeometryModel`](@ref) from a hierarchy of volumes. Boolean shapes are converted to
CSG programs over primitive shapes, and a BVH is built for the daughters of volumes with more than
`bvh_threshold` daughters.
"""
function close_geometry(world::Volume{T}; bvh_threshold::Integer=4) where T
    b = _ModelBuilder{T}(bvh_threshold)
    _add_volume!(b, world, 0)
    shapes = Tuple(b.arrays)
    D = isempty(b.depth) ? 0 : maximum(values(b.depth))
    GeometryModel{T,D}(shapes, b.volumes, b.placements, b.materials, b.bvhnodes, b.bvhindices,
                       b.cones, b.zoffsets, b.csg_components, b.csg_program)
end
