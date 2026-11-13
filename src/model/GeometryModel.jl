#---Flat geometry model-------------------------------------------------------------------------------
#   A "closed" geometry (see `close_geometry`) where the hierarchy of volumes is stored in flat arrays
#   of isbits structures, referenced by indices. The same model (and the same navigation code) can be
#   used on the CPU or, after `adapt`ing the arrays, on a GPU.
#
#   Shapes are stored in one array per concrete type (the `shapes` tuple) and referenced by a
#   `ShapeRef(kind, index)`. Calling a shape function goes through `dispatch`, a generated chain of
#   `if` statements over the kinds, so there is no dynamic dispatch. Polycones keep their sections in
#   shared arrays. Boolean shapes (of any nesting) are stored as flat CSG expressions over primitive
#   shapes (see CSG.jl), so that the generated code does not grow with the nesting of booleans: some
#   GPU compilers (e.g. Metal) inline all function calls.
#---------------------------------------------------------------------------------------------------

"Reference to a shape of a `GeometryModel`: kind (index in the tuple of shape arrays) and index."
struct ShapeRef
    kind::Int32
    index::Int32
end

struct FlatVolume
    shape::ShapeRef
    material::Int32
    first_daughter::Int32        # index of the first daughter in `placements` minus one
    ndaughters::Int32
    bvh::Int32                   # first node of the BVH of the daughters (0 if none)
end

struct FlatPlacement{T<:AbstractFloat}
    volume::Int32                          # placed (logical) volume
    transformation::Transformation3D{T}    # mother frame -> daughter frame
    aabb::AABB{T}                          # bounding box of the daughter in the mother frame
end

struct FlatMaterial{T<:AbstractFloat}
    density::T
    temperature::T
    Amass::T
end

"Node of a BVH stored in depth-first order with escape indices (stackless traversal)."
struct BVHNode{T<:AbstractFloat}
    aabb::AABB{T}
    first::Int32     # leaf: first entry in `bvhindices` (0 for internal nodes)
    count::Int32     # leaf: number of entries
    skip::Int32      # next node when the subtree is not visited (0 = end of traversal)
end

#---Flat representations of shapes that have components----------------------------------------------
struct FlatPolycone{T<:AbstractFloat}
    first::Int32     # offset of the sections in `cones` / `zoffsets`
    count::Int32
    ϕ₀::T
    Δϕ::T
end

"Component of a CSG expression: a primitive shape placed in the frame of the boolean shape."
struct CSGComponent{T<:AbstractFloat}
    ref::ShapeRef                          # primitive shape
    transformation::Transformation3D{T}    # CSG frame -> primitive frame
    aabb::AABB{T}                          # bounding box in the CSG frame (with tolerance)
end

# Operations of the CSG programs (positive values push the inside-ness of a component)
const CSG_UNION        = Int32(-1)
const CSG_SUBTRACTION  = Int32(-2)
const CSG_INTERSECTION = Int32(-3)

"Boolean shape as a postfix program over `ncomp` components (stored in shared arrays of the model)."
struct FlatCSG{T<:AbstractFloat}
    first_comp::Int32    # offset in `csg_components`
    ncomp::Int32
    first_op::Int32      # offset in `csg_program`
    nops::Int32
    aabb::AABB{T}        # extent of the boolean shape
end

"Shape that is never inside and never hit (replaces unsupported shapes)."
struct FlatNoShape{T<:AbstractFloat} <: AbstractShape{T} end
inside(::FlatNoShape{T}, ::Point3{T}) where T = kOutside
distance_to_in(::FlatNoShape{T}, ::Point3{T}, ::Vector3{T}) where T = T(Inf)
distance_to_out(::FlatNoShape{T}, ::Point3{T}, ::Vector3{T}) where T = T(-1)
safety_to_in(::FlatNoShape{T}, ::Point3{T}) where T = T(Inf)
safety_to_out(::FlatNoShape{T}, ::Point3{T}) where T = T(-1)
normal(::FlatNoShape{T}, ::Point3{T}) where T = zero(Vector3{T})
extent(::FlatNoShape{T}) where T = (zero(Point3{T}), zero(Point3{T}))

#---The model-----------------------------------------------------------------------------------------
"""
    GeometryModel

Flat representation of a geometry, built with [`close_geometry`](@ref). `D` is 1 if the model has
boolean (CSG) shapes, 0 otherwise. Use `Adapt.adapt(backend, model)` to move it to a GPU.
"""
struct GeometryModel{T<:AbstractFloat, D, S<:Tuple, VV, VP, VM, VN, VI, VC, VZ, VK, VO}
    shapes::S                     # tuple of arrays, one per kind of shape
    volumes::VV                   # FlatVolume
    placements::VP                # FlatPlacement{T}
    materials::VM                 # FlatMaterial{T}
    bvhnodes::VN                  # BVHNode{T}
    bvhindices::VI                # Int32 (placement indices)
    cones::VC                     # Cone{T}: sections of polycones
    zoffsets::VZ                  # T: z offsets of the sections of polycones
    csg_components::VK            # CSGComponent{T}: components of boolean shapes
    csg_program::VO               # Int32: programs of boolean shapes
end

function GeometryModel{T,D}(shapes::S, volumes::VV, placements::VP, materials::VM, bvhnodes::VN, bvhindices::VI,
                            cones::VC, zoffsets::VZ, csg_components::VK, csg_program::VO) where {T,D,S,VV,VP,VM,VN,VI,VC,VZ,VK,VO}
    GeometryModel{T,D,S,VV,VP,VM,VN,VI,VC,VZ,VK,VO}(shapes, volumes, placements, materials, bvhnodes, bvhindices,
                                                    cones, zoffsets, csg_components, csg_program)
end

function Adapt.adapt_structure(to, m::GeometryModel{T,D}) where {T,D}
    GeometryModel{T,D}(map(a -> adapt(to, a), m.shapes), adapt(to, m.volumes), adapt(to, m.placements),
                       adapt(to, m.materials), adapt(to, m.bvhnodes), adapt(to, m.bvhindices),
                       adapt(to, m.cones), adapt(to, m.zoffsets), adapt(to, m.csg_components), adapt(to, m.csg_program))
end

function Base.show(io::IO, m::GeometryModel{T,D}) where {T,D}
    print(io, "GeometryModel{$T}(volumes=", length(m.volumes), ", placements=", length(m.placements),
          ", materials=", length(m.materials), ", shape kinds=", length(m.shapes),
          ", shapes=", sum(length, m.shapes), ", bvh nodes=", length(m.bvhnodes), ", csg components=", length(m.csg_components), ")")
end

#---Materialization of flat shapes into shapes with the generic algorithms----------------------------
@inline _materialize(model, shape, ::Val) = shape
@inline function _materialize(model::GeometryModel{T}, s::FlatPolycone{T}, ::Val) where T
    r = s.first+1:s.first+s.count
    Polycone{T}(view(model.cones, r), view(model.zoffsets, r), s.ϕ₀, s.Δϕ)
end
@inline function _materialize(model::M, s::FlatCSG{T}, ::Val{D}) where {M,T,D}
    D == 0 && return FlatNoShape{T}()     # CSG components are primitives (never happens)
    CSGShape{T,M}(model, s)
end

"""
    dispatch(f, model, ref, Val(D), args...)

Call `f(shape, args...)` for the shape referenced by `ref`. `D = 1` for the shapes of volumes, which
can be boolean (CSG) shapes, `D = 0` for the primitive components of CSG shapes.
"""
@generated function _dispatch(f::F, model::GeometryModel{T,DM,S}, ref::ShapeRef, ::Val{D}, args...) where {F,T,DM,S,D}
    N = length(S.parameters)
    body = :(f(_materialize(model, @inbounds(model.shapes[$N][ref.index]), Val(D)), args...))
    for i in N-1:-1:1
        body = :(ref.kind == $i ? f(_materialize(model, @inbounds(model.shapes[$i][ref.index]), Val(D)), args...) : $body)
    end
    return body
end
# Inlined: twice faster on the CPU; GPU compilers such as Metal's inline all calls anyway (the CSG
# representation of booleans keeps the size of the generated code linear in the number of kinds).
@inline dispatch(f::F, model::GeometryModel, ref::ShapeRef, ::Val{D}, args...) where {F,D} = _dispatch(f, model, ref, Val(D), args...)

#---Shape functions on the shape of a volume of the model----------------------------------------------
@inline shape_call(f, model::GeometryModel{T,D}, ref::ShapeRef, args...) where {T,D} = dispatch(f, model, ref, Val(D), args...)
