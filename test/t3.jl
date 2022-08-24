using Revise
using AbstractTrees
using BenchmarkTools

struct BVH
    name::String
    children::Union{Tuple{}, Tuple{BVH, BVH}}
    indices::Vector{Int64}
end

BVH(n::String, l::BVH, r::BVH) = BVH(n, (l,r), Int64[])   # construct a node
BVH(n::String, inds::Vector{Int64}) = BVH(n, (), inds)    # construct a leaf

Base.show(io::IO, b::BVH) = print(io, isempty(b.children) ? "Leaf-" : "Node-", b.name)
AbstractTrees.children(bvh::BVH) = bvh.children

function _indices!(ind::Vector{Int64}, b::BVH)
    isempty(b.children) && return append!(ind, b.indices)
    _indices!(ind, b.children[1])
    _indices!(ind, b.children[2])
    return ind
end
indices(b::BVH) = _indices!(Int[], b)
indices!(ind::Vector{Int64}, b::BVH) = return _indices!(empty!(ind), b)


bvh = BVH("root", BVH("1",BVH("2", [1,2,3]),BVH("3", BVH("4", [4,5,6]), BVH("5",[7,8,9]))),BVH("6",[10,11,12]))

print_tree(bvh)
collect(indices(bvh))

@btime indices($bvh)

ind = zeros(Int64,10000)
@btime indices!($ind, $bvh);

function f(b::BVH)
    sum = 0
    for i in indices!(ind, b)
        sum += i
    end
    return sum
end
@btime f($bvh)



#------
using AbstractTrees
using BenchmarkTools

struct BVH
    name::String
    children::Union{Nothing, Tuple{BVH, BVH}}
    indices::Union{Nothing, Vector{Int64}}
end

BVH(n::String, l::BVH, r::BVH) = BVH(n, (l,r), nothing)        # construct a node
BVH(n::String, inds::Vector{Int64}) = BVH(n, nothing, inds)    # construct a leaf

Base.show(io::IO, b::BVH) = print(io, isnothing(b.children) ? "Leaf-" : "Node-", b.name)
AbstractTrees.children(bvh::BVH) = isnothing(bvh.children) ? [] : bvh.children

function indices(b::BVH)
    isnothing(b.children) && return (i for i in b.indices)
    isnothing(b.indices) && return Iterators.flatten(indices(c) for c in b.children)
end

bvh = BVH("root", BVH("1",BVH("2", [1,2,3]),BVH("3", BVH("4", [4,5,6]), BVH("5",[7,8,9]))),BVH("6",[10,11,12]))

@btime indices($bvh)


function f(b::BVH)
    sum = 0
    for i in indices(b)
        sum += i
    end
    return sum
end



@btime f($bvh)
