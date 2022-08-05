#---Boxed volumes used to build the Bounding Volume Hierarchy---------------------------------------
#   BoxedVolumes wraps around a set of PlacedVolumes indices, and the AABBs (axis aligned bounding 
#   boxes) of the corresponding daughters, identified by min/max coordinates

struct BoxedVolumes{T<:AbstractFloat}
    min::Matrix{T}
    max::Matrix{T}
    center::Matrix{T}
    indices::Vector{Int}
end

#---BoxedVoumes constructor using the lists f daughters----------------------------------------------
function BoxedVolumes(pvols::Vector{PlacedVolume{T}}) where T
    min_ = [aabb.min[i] for aabb in (AABB(pvol) for pvol in  pvols), i in 1:3]'
    max_ = [aabb.max[i] for aabb in (AABB(pvol) for pvol in  pvols), i in 1:3]'
    center_ = 0.5 * (min_ + max_)
    BoxedVolumes{T}(min_, max_, center_, 1:length(pvols))
end

Base.length(bv::BoxedVolumes) = length(bv.indices)

function Base.sort!(bv::BoxedVolumes, i::Int)
    perm          = sortperm(bv.center[i, :])
    bv.min[:]     = bv.min[:, perm]
    bv.max[:]     = bv.max[:, perm]
    bv.center[:]  = bv.center[:, perm]
    bv.indices[:] = bv.indices[perm]
    nothing
end

function Base.reverse!(bv::BoxedVolumes)
    bv.min[:] = bv.min[:, end:-1:1]
    bv.max[:] = bv.max[:, end:-1:1]
    bv.center[:] = bv.center[:, end:-1:1]
    bv.indices[:] = bv.indices[end:-1:1]
    nothing
end

function set_area!(area::Vector{T}, bv::BoxedVolumes{T})  where T
    tmpmin, tmpmax = bv.min[:,1], bv.max[:,1]
    area[1] = _area(tmpmax .- tmpmin)
    for i in 2:(length(bv)-1)
        tmpmin[:] = min.(tmpmin, bv.min[:,i])
        tmpmax[:] = max.(tmpmax, bv.max[:,i])
        area[i] = _area(tmpmax .- tmpmin)
    end
    nothing
end

function optimal_split(bv::BoxedVolumes{T}, axis::Int, buffer1::Vector{T}, buffer2::Vector{T}) where T
    N = length(bv)
    @assert N-1 == length(buffer1) == length(buffer2)
    sort!(bv, axis)
    set_area!(buffer1, bv)
    reverse!(bv)
    set_area!(buffer2, bv)
    return minimum((cost=i*a1_ + (N-i)*a2_, index=i, sortaxis=axis)
                   for (i, (a1_, a2_)) in enumerate(zip(buffer1,reverse(buffer2))))
end

function optimal_split(bv::BoxedVolumes{T}) where T
    N = length(bv)
    # create buffers
    buffer1 = Vector{T}(undef, N-1)
    buffer2 = Vector{T}(undef, N-1)
    min(optimal_split(bv, 1, buffer1, buffer2),
        optimal_split(bv, 2, buffer1, buffer2),
        optimal_split(bv, 3, buffer1, buffer2))
end

function subdivide(bv::BoxedVolumes, optsplit)
    i, axis = optsplit.index, optsplit.sortaxis
    sort!(bv, axis)
    (BoxedVolumes(bv.min[:,1:i], bv.max[:,1:i], bv.center[:,1:i], bv.indices[1:i]),
     BoxedVolumes(bv.min[:,i+1:end],bv.max[:,i+1:end], bv.center[:,i+1:end], bv.indices[i+1:end]))
end

function AABB(bv::BoxedVolumes{T}) where T
    AABB{T}(minimum(bv.min, dims=2) .- eps(T), maximum(bv.max, dims=2) .+ eps(T))
end

#---BVH---------------------------------------------------------------------------------------------
struct BVHParam{T<:AbstractFloat}
    SAHfrac::T
    BVHParam{T}(SAHfrac=T(8)) where T = new(SAHfrac)
end

struct BVH{T<:AbstractFloat,N,C}
    aabb::AABB{T}
    children::NTuple{N, C} # either N=2 and C=BVH, or N=1 and C=PVolIndices
    BVH(bv::AABB{T}, children...) where T = new{T, length(children), eltype(children)}(bv, children)
end
children(n::BVH) = n.children

struct PVolIndices
    inds::Vector{Int64}
end
children(n::PVolIndices) = ()

function Base.show(io::IO, bvh::BVH{T, N, C}) where {T,N,C}
    print(io, "BVH{$T}",(aabb=bvh.aabb, leaves=length(collect(pvolindices(bvh)))))
end


#--- _Branch and _Leaf are only used to construct the BVH recursively
struct _Branch{T} data::T; end
struct _Leaf{T} data::T; end

_children(::_Leaf, param::BVHParam) = ()
function _children(n::_Branch, bp::BVHParam)
    bc = n.data
    N, A = length(bc), area(AABB(bc))
    if N>1
        optsplit = optimal_split(bc)
        # optsplit.cost / A approximates the expected number of checks when partitioning the cells
        if  N - optsplit.cost/A > bp.SAHfrac
            bc1, bc2 = subdivide(bc, optsplit)
            return (_Branch(bc1), _Branch(bc2))
        end
    end
    return (_Leaf(bc),)
end
_buildBVH(n::_Leaf, ::BVHParam) = PVolIndices(n.data.indices)
_buildBVH(n::_Branch, bp::BVHParam) = BVH(AABB(n.data), [_buildBVH(c, bp) for c in _children(n, bp)]...)

buildBVH(pvols::Vector{PlacedVolume{T}}, bp::BVHParam{T}) where T = _buildBVH(_Branch(BoxedVolumes(pvols)), bp)
buildBVH(pvols::Vector{PlacedVolume{T}}) where T = buildBVH(pvols, BVHParam{T}())

resetBVH(bvh::BVH, vertices, faces) = nothing

#---Returns an iterator over all subsets of PlacedVolume indices such that `f(x::AABB) == true` for all parent axis aligned bounding boxes `x`.

pvolindices(f, n::BVH{T, N, C}) where {T<:AbstractFloat, N, C<:PVolIndices} = (c.inds for c in n.children)
pvolindices(f, n::BVH{T, N, C}) where {T<:AbstractFloat, N, C<:BVH} = Iterators.flatten(pvolindices(f, c) for c in children(n) if isa(c, PVolIndices) || f(c.aabb))
#pvolindices(f, n::BVH{T, N, C}) where {T<:AbstractFloat, N, C<:BVH} = Iterators.flatten(i for c in children(n) if isa(c, PVolIndices) || f(c.aabb) for i in c)
pvolindices(head::BVH{T, N, C}) where {T<:AbstractFloat, N, C<:BVH} = pvolindices(x->true, head)
pvolindices(head::BVH{T, N, C}) where {T<:AbstractFloat, N, C<:PVolIndices} = (c.inds for c in head.children)

function _pushPvolIndices!(ind::Vector{Int64}, f::Function, b::BVH{T, N, C}) where {T<:AbstractFloat, N, C<:BVH}
    for c in children(b)
        if f(c.aabb)
            _pushPvolIndices!(ind, f, c)
        end
    end
end
function _pushPvolIndices!(ind::Vector{Int64}, ::Function, b::BVH{T, N, C}) where {T<:AbstractFloat, N, C<:PVolIndices}
    append!(ind, b.children[1].inds)
end
function pushPvolIndices!(ind::Vector{Int64}, f::Function, b::BVH{T, N, C}) where {T<:AbstractFloat, N, C}
    empty!(ind)
    _pushPvolIndices!(ind, f, b)
end




#=

function insideBVH(point::Point3{T}, vol::Volume{T}, bvh::BVH{T})
    pvols = vol.daughters
    for inds in pvolindices( x-> inside(point,x), bvh)
        for i in inds
            @show i
            Geom4hep.contains(pvols[i], point) && return true
        end
    end
    return false
end

function insideNaive(point::Point3{T}, vol::Volume{T})
    pvols = vol.daughters
    for pvol in pvols
        Geom4hep.contains(pvol, point) && return true
    end
    return false
end

=#