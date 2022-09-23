#---Boolean Types----------------------------------------------------------------------------------
abstract type AbstractBoolean{T, SL, SR} <: AbstractShape{T} end

struct BooleanUnion{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape} <: AbstractBoolean{T, SL, SR}
    left::SL        # the mother (or left) volume A in unplaced form
    right::SR       # (or right) volume B in placed form, acting on A with a boolean operation
    transformation::Transformation3D{T} # placement of "right" with respect of "left"
end
struct BooleanSubtraction{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape} <: AbstractShape{T}
    left::SL        # the mother (or left) volume A in unplaced form
    right::SR       # (or right) volume B in placed form, acting on A with a boolean operation
    transformation::Transformation3D{T} # placement of "right" with respect of "left"
end
struct BooleanIntersection{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape} <: AbstractShape{T}
    left::SL        # the mother (or left) volume A in unplaced form
    right::SR       # (or right) volume B in placed form, acting on A with a boolean operation
    transformation::Transformation3D{T} # placement of "right" with respect of "left"
end


#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, shape::BooleanUnion{T, SL, SR}) where {T,SL,SR}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanUnion",(left=left, right=right, placement=placement))
end
function Base.show(io::IO, shape::BooleanSubtraction{T, SL, SR}) where {T,SL,SR}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanSubtraction",(left=left, right=right, placement=placement))
end
function Base.show(io::IO, shape::BooleanIntersection{T, SL, SR}) where {T,SL,SR}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanIntersection",(left=left, right=right, placement=placement))
end


function GeometryBasics.mesh(shape::AbstractBoolean{T, SL, SR}) where {T,SL,SR}
    (; left, right, transformation) = shape
    rmesh = mesh(right)
    merge([mesh(left), Mesh(map(c -> Point3{T}(c * transformation), coordinates(rmesh)), faces(rmesh))])
end

#---Basic functions---------------------------------------------------------------------------------
function extent(shape::BooleanUnion{T, SL, SR})::Tuple{Point3{T},Point3{T}} where {T,SL,SR}
    (; left, right, transformation) = shape
    minLeft, maxLeft = extent(left)
    minRight, maxRight = extent(right) .* Ref(transformation)
    (min.(minLeft, minRight), max.(maxLeft, maxRight))
end
function extent(shape::BooleanSubtraction{T, SL, SR})::Tuple{Point3{T},Point3{T}} where {T,SL,SR}
    extent(shape.left)
end
function extent(shape::BooleanIntersection{T, SL, SR})::Tuple{Point3{T},Point3{T}} where {T,SL,SR}
    (; left, right, transformation) = shape
    minLeft, maxLeft = extent(left)
    minRight, maxRight = extent(right) .* Ref(transformation)
    (max.(minLeft, minRight), min.(maxLeft, maxRight))
end

 
function safetyToOut(shape::AbstractBoolean)
    return 0
end
function safetyToIn(shape::AbstractBoolean)
    return 0
end
