#---Boolean-------------------------------------------------------------
struct Boolean{T<:AbstractFloat, SL<:AbstractShape, SR<:AbstractShape, OP} <: AbstractShape{T}
    left::SL        # the mother (or left) volume A in unplaced form
    right::SR       # (or right) volume B in placed form, acting on A with a boolean operation
    transformation::Transformation3D{T} # placement of "right" with respect of "left"
end

#---Constructor------------------------------------------------------------------------------------
function Boolean(op::Symbol, left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    @assert op in (:Union, :Subtraction, :Intersection) "Boolean suppoted operations are: :Union, :Subtraction, :Intersection"
    Boolean{T,typeof(left),typeof(right),op}(left, right, place)
end
#function BooleanUnion(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
#    Boolean{T,typeof(left),typeof(right),:Union}(left, right, place)
#end
function Subtraction(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    Boolean{T,typeof(left),typeof(right),:Subtraction}(left, right, place)
end
function Intersection(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    Boolean{T,typeof(left),typeof(right),:Intersection}(left, right, place)
end

#---Utilities---------------------------------------------------------------------------------------

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, shape::Boolean{T, SL, SR, OP}) where {T,SL,SR,OP}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "Boolean{$OP}",(left=left, right=right, placement=placement))
end

function GeometryBasics.coordinates(shape::Boolean{T, SL, SR, OP}, facets=36) where {T,SL,SR,OP}
    (; left, right, transformation) = shape
    coors = (coordinates(left, facets), coordinates(right, facets))
    transform = (one(Transformation3D{T}), transformation)
    return (c * transform[i] for i in 1:2 for c in coors[i])
end

function GeometryBasics.faces(shape::Boolean{T, SL, SR, OP}, facets=36) where {T,SL,SR,OP}
    (; left, right) = shape
    offset = length(coordinates(left, facets))
    append!(faces(left, facets), [t .+ offset for t in faces(right, facets)])
end

#---Basic functions---------------------------------------------------------------------------------
function extent(shape::Boolean{T, SL, SR, OP})::Tuple{Point3{T},Point3{T}} where {T,SL,SR,OP}
    (; left, right, transformation) = shape
    minLeft, maxLeft = extent(left)
    minRight, maxRight = extent(right) .* Ref(transformation)
    if OP == :Union
        (min.(minLeft, minRight), max.(maxLeft, maxRight))
    elseif  OP == :Subtraction
        (minLeft, maxLeft)
    elseif  OP == :Intersection
        (max.(minLeft, minRight), min.(maxLeft, maxRight))
    end
end
 
function inside(shape::Boolean{T, SL, SR, :Union}, point::Point3{T}) where {T,SL,SR}
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

function inside(shape::Boolean{T, SL, SR, :Intersection}, point::Point3{T}) where {T,SL,SR}
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

function inside(shape::Boolean{T, SL, SR, :Subtraction}, point::Point3{T}) where {T,SL,SR}
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

function safetyToOut(shape::Boolean{T, SL, SR, OP}, point::Point3{T})::T where {T,SL,SR, OP}
    return 0
end
function safetyToIn(shape::Boolean{T, SL, SR, OP}, point::Point3{T})::T where {T,SL,SR, OP}
    return 0
end

function distanceToOut(shape::Boolean{T, SL, SR, :Union}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    dist = T(0)
    positionA = inside(left, point)
    if positionA != kOutside  # point inside A
        while(true)
            distA = distanceToOut(left, point, dir)
            dist += (distA > 0 && distA < Inf) ? distA : 0 + kPushTolerance(T)
            npoint = point + dist * dir
            lpoint = transformation * npoint
            if inside(right, lpoint) != kOutside # B could be overlapping with A -- and/or connecting A to another part of A
                ldir = transformation * dir
                distB = distanceToOut(right, lpoint, ldir)
                dist += (distB > 0 && distB < Inf) ? distB : 0 + kPushTolerance(T)
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
                dist += (distB > 0 && distB < Inf) ? distB : 0 + kPushTolerance(T)
                npoint = point + dist * dir
                if inside(left, npoint) != kOutside # A could be overlapping with B -- and/or connecting B to another part of B
                    ldir = transformation * dir
                    distA = distanceToOut(left, npoint, dir)
                    dist += (distA > 0 && distA < Inf) ? distA : 0 + kPushTolerance(T)
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

function distanceToOut(shape::Boolean{T, SL, SR, :Intersection}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    distA = distanceToOut(left, point, dir)
    distB = distanceToOut(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distanceToOut(shape::Boolean{T, SL, SR, :Subtraction}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    distA = distanceToOut(left, point, dir)
    distB = distanceToIn(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distanceToIn(shape::Boolean{T, SL, SR, :Union}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
    (; left, right, transformation) = shape
    distA = distanceToIn(left, point, dir)
    distB = distanceToIn(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distanceToIn(shape::Boolean{T, SL, SR, :Intersection}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
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

function distanceToIn(shape::Boolean{T, SL, SR, :Subtraction}, point::Point3{T}, dir::Vector3{T})::T where {T,SL,SR}
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
            npoint = point + dist * dir
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
        dist += (d1 >= 0 && d1 < Inf) ? d1 + kPushTolerance(T) : 0
        point + dist * dir
        lpoint = transformation * npoint
        inRight = true
    end
end
