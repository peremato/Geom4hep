@compactify begin
    @abstract struct AbstractBoolean{T, SL, SR} <: AbstractShape{T}
        left::SL = NoShape()
        right::SR = NoShape()
        transformation::Transformation3D{T} = Transformation3D{T}(0, 0, 0)
    end
    struct BooleanUnion{T, SL, SR} <: AbstractBoolean{T, SL, SR} end
    struct BooleanSubtraction{T, SL, SR} <: AbstractBoolean{T, SL, SR} end
    struct BooleanIntersection{T, SL, SR} <: AbstractBoolean{T, SL, SR} end
end
#---Constructor------------------------------------------------------------------------------------
function BooleanUnion(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanUnion{T,typeof(left),typeof(right)}(; left, right, transformation = place)
end
function BooleanSubtraction(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanSubtraction{T,typeof(left),typeof(right)}(; left, right, transformation = place)
end
function BooleanIntersection(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanIntersection{T,typeof(left),typeof(right)}(; left, right, transformation = place)
end

#---Boolean Types----------------------------------------------------------------------------------
# struct BooleanUnion{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape} <: AbstractShape{T}
#     left::SL        # the mother (or left) volume A in unplaced form
#     right::SR       # (or right) volume B in placed form, acting on A with a boolean operation
#     transformation::Transformation3D{T} # placement of "right" with respect of "left"
# end
# struct BooleanSubtraction{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape} <: AbstractShape{T}
#     left::SL        # the mother (or left) volume A in unplaced form
#     right::SR       # (or right) volume B in placed form, acting on A with a boolean operation
#     transformation::Transformation3D{T} # placement of "right" with respect of "left"
# end
# struct BooleanIntersection{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape} <: AbstractShape{T}
#     left::SL        # the mother (or left) volume A in unplaced form
#     right::SR       # (or right) volume B in placed form, acting on A with a boolean operation
#     transformation::Transformation3D{T} # placement of "right" with respect of "left"
# end

#---Utilities---------------------------------------------------------------------------------------

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


function GeometryBasics.mesh(shape::BooleanUnion{T, SL, SR}) where {T,SL,SR}
    (; left, right, transformation) = shape
    rmesh = mesh(right)
    merge([mesh(left), Mesh(map(c -> Point3{T}(c * transformation), coordinates(rmesh)), faces(rmesh))])
end
function GeometryBasics.mesh(shape::BooleanSubtraction{T, SL, SR}) where {T,SL,SR}
    (; left, right, transformation) = shape
    rmesh = mesh(right)
    merge([mesh(left), Mesh(map(c -> Point3{T}(c * transformation), coordinates(rmesh)), faces(rmesh))])
end
function GeometryBasics.mesh(shape::BooleanIntersection{T, SL, SR}) where {T,SL,SR}
    (; left, right, transformation) = shape
    rmesh = mesh(right)
    merge([mesh(left), Mesh(map(c -> Point3{T}(c * transformation), coordinates(rmesh)), faces(rmesh))])
end

#---Basic functions---------------------------------------------------------------------------------
function extent(shape::AbstractBoolean)
    (; left, right, transformation) = shape
    @compactified shape::AbstractBoolean begin
        BooleanUnion => begin
            minLeft, maxLeft = extent(left)
            minRight, maxRight = extent(right) .* Ref(transformation)
            (min.(minLeft, minRight), max.(maxLeft, maxRight))
        end
        BooleanIntersection => begin
            minLeft, maxLeft = extent(left)
            minRight, maxRight = extent(right) .* Ref(transformation)
            (max.(minLeft, minRight), min.(maxLeft, maxRight))
        end
        BooleanSubtraction => extent(shape.left)
    end
end

 
function safetyToOut(shape::AbstractBoolean)
    return 0
end
function safetyToIn(shape::AbstractBoolean)
    return 0
end
