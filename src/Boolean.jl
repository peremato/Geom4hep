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

function distanceToOut(shape::BooleanUnion{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    dist = T(0)
    positionA = inside(left, point)
    if positionA != kOutside  # point inside A
        while(true)
            distA = distanceToOut(left, point, dir)
            dist += (distA > 0 && distA < Inf) ? distA : 0
            dist += kPushTolerance(T)
            npoint = point + dist * dir
            lpoint = transformation * npoint
            if inside(right, lpoint) != kOutside # B could be overlapping with A -- and/or connecting A to another part of A
                ldir = transformation * dir
                distB = distanceToOut(right, lpoint, ldir)
                dist += (distB > 0 && distB < Inf) ? distB : 0
                dist += kPushTolerance(T)
                npoint = point + dist * dir
                if inside(left, npoint) == kOutside
                    break
                end
            else
                break
            end
        end
        return dist - kPushTolerance(T)
    else
        lpoint = transformation * point
        positionB = inside(right, lpoint)
        if positionB != kOutside  # point inside B
            ldir = transformation * dir
            while(true)
                distB = distanceToOut(right, lpoint, ldir)
                dist += (distB > 0 && distB < Inf) ? distB : 0
                dist += kPushTolerance(T) # Give a push
                npoint = point + dist * dir
                if inside(left, npoint) != kOutside # A could be overlapping with B -- and/or connecting B to another part of B
                    ldir = transformation * dir
                    distA = distanceToOut(left, npoint, dir)
                    dist += (distA > 0 && distA < Inf) ? distA : 0
                    dist += kPushTolerance(T)
                    npoint = point + dist * dir
                    lpoint = transformation * npoint
                    if inside(right, lpoint) == kOutside
                        break
                    end
                else
                    break
                end
            end
            return dist - kPushTolerance(T)
        else
            return T(-1)
        end
    end 
end

function distanceToOut(shape::BooleanIntersection{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    distA = distanceToOut(left, point, dir)
    distB = distanceToOut(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distanceToOut(shape::BooleanSubtraction{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    distA = distanceToOut(left, point, dir)
    distB = distanceToIn(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distanceToIn(shape::BooleanUnion{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    distA = distanceToIn(left, point, dir)
    distB = distanceToIn(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distanceToIn(shape::BooleanIntersection{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape

    positionA = inside(left, point)
    lpoint = transformation * point
    positionB = inside(right, lpoint)

    inLeft = positionA == kInside
    inRight = positionB == kInside

    inLeft && inRight && return T(-1)

    dist = T(0)
    npoint = point
    ldir = transformation * dir
    #  main loop
    while true
        d1 = d2 = 0
        if !inLeft
            d1 = distanceToIn(left, npoint, dir)
            d1 = max(d1, kTolerance(T))
            d1 == T(Inf) && return T(Inf)
        end
        if !inRight
            d2 = distanceToIn(right, lpoint, ldir)
            d2 = max(d2, kTolerance(T))
            d2 == T(Inf) && return T(Inf)
        end
        if d1 > d2
            # propagate to left shape
            dist += d1
            inleft = true
            npoint += d1 * dir
            lpoint = transformation * npoint
            # check if propagated point is inside right shape, check is done with a little push
            inRight = inside(right, lpoint + kTolerance(T) * ldir) == kInside
            inRight && return dist
            # here inleft=true, inright=false
        else
            # propagate to right shape
            dist += d2
            inright = true
            # check if propagated point is inside left shape, check is done with a little push
            npoint += d2 * dir
            lpoint = transformation * npoint
            inLeft = inside(left, npoint + kTolerance(T) * dir) == kInside
            inLeft && return dist
        end
    end
    return dist
end

function distanceToIn(shape::BooleanSubtraction{T, SL, SR}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape

    lpoint = transformation * point
    ldir = transformation * dir
    positionB = inside(left, lpoint)
    inRight = positionB == kInside

    npoint = point
    dist = T(0)

    while true
        if inRight
            # propagate to outside of '- / RightShape'
            d1 = distanceToOut(right, lpoint, ldir)
            dist += (d1 >= 0 && d1 < Inf) ? d1 + kPushTolerance(T) : 0
            npoint = point + (dist + kPushTolerance(T)) * dir
            lpoint = transformation * npoint
            # now master outside 'B'; check if inside 'A'
            inside(left, npoint) == kInside && distanceToOut(left, npoint, dir) > kPushTolerance(T) && return dist
        end

        # if outside of both we do a max operation master outside '-' and outside '+' ;  find distances to both
        d2 = distanceToIn(left, npoint, dir)
        d2 = max(d2, 0)
        d2 == T(Inf) && return T(Inf)

        d1 = distanceToIn(right, lpoint, ldir)
        if d2 < d1 - kTolerance(T)
          dist += d2 + kPushTolerance(T)
          return dist
        end

        #   propagate to '-'
        dist += (d1 >= 0 && d1 < Inf) ? d1 : 0
        npoint = point + (dist + kPushTolerance(T)) * dir
        lpoint = transformation * npoint
        inRight = true
    end
end
