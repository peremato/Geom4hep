#---Boolean Types----------------------------------------------------------------------------------
struct BooleanUnion{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape} <: AbstractShape{T}
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

#---Constructor------------------------------------------------------------------------------------
function BooleanUnion(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanUnion{T,typeof(left),typeof(right)}(left, right, place)
end
function BooleanSubtraction(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanSubtraction{T,typeof(left),typeof(right)}(left, right, place)
end
function BooleanIntersection(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanIntersection{T,typeof(left),typeof(right)}(left, right, place)
end

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

 
function inside(shape::BooleanUnion{T, SL, SR}, point::Point3{T})::Int64 where {T,SL,SR}
    (; left, right, transformation) = shape
    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kInside && return kInside

    positionB = inside(right, lpoint)
    positionB == kInside && return kInside

    if positionA == kSurface && positionB == kSurface
        normalA = normal(left, point)
        normalB = normal(right, lpoint) * transformation
        if dot(normalA, normalB) < 0
            return kInside    # touching solids -)(-
        else 
            return kSurface   # overlapping solids =))
        end
    elseif positionA == kSurface || positionB == kSurface
        return kSurface
    else
        return kOutside
    end
end

function inside(shape::BooleanIntersection{T, SL, SR}, point::Point3{T})::Int64  where {T,SL,SR}
    (; left, right, transformation) = shape
    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kOutside && return kOutside

    positionB = inside(right, lpoint)
    if positionA == kInside && positionB == kInside
        return kInside
    elseif  positionA == kInside && positionB == kSurface ||
            positionA == kSurface && positionB == kInside ||
            positionA == kSurface && positionB == kSurface
        return kSurface
    else
        return kOutside
    end 
end

function inside(shape::BooleanSubtraction{T, SL, SR}, point::Point3{T})::Int64  where {T,SL,SR}
    (; left, right, transformation) = shape
    lpoint = transformation * point

    positionA = inside(left, point)
    positionA == kOutside && return kOutside

    positionB = inside(right, lpoint)

    if positionA == kInside && positionB == kOutside
        return kInside;
    elseif positionA == kInside && positionB == kSurface ||
           positionB == kOutside && positionA == kSurface
        return kSurface
    else
        return kOutside
    end
end

function safetyToOut(shape::BooleanUnion{T, SL, SR}, point::Point3{T})::T where {T,SL,SR}
    return 0
end
function safetyToOut(shape::BooleanSubtraction{T, SL, SR}, point::Point3{T})::T where {T,SL,SR}
    return 0
end
function safetyToOut(shape::BooleanIntersection{T, SL, SR}, point::Point3{T})::T where {T,SL,SR}
    return 0
end

function safetyToIn(shape::BooleanUnion{T, SL, SR}, point::Point3{T})::T where {T,SL,SR}
    return 0
end
function safetyToIn(shape::BooleanSubtraction{T, SL, SR}, point::Point3{T})::T where {T,SL,SR}
    return 0
end
function safetyToIn(shape::BooleanIntersection{T, SL, SR}, point::Point3{T})::T where {T,SL,SR}
    return 0
end
