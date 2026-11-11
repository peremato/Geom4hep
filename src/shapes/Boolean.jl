#---Basic Shape--------------------------------------------------------------------------------------
const BaseShape{T<:AbstractFloat} = Union{NoShape{T},
                                          Box{T},
                                          Trd{T},
                                          Trap{T},
                                          Tube{T},
                                          Cone{T},
                                          Polycone{T},
                                          CutTube{T}}

#---Boolean Shapes----------------------------------------------------------------------------------
struct BooleanShape{T<:AbstractFloat, OP} <: AbstractShape{T}
    left::Union{BaseShape{T}, BooleanShape{T}}   # the mother (or left) volume A in unplaced form
    right::Union{BaseShape{T}, BooleanShape{T}}  # (or right) volume B in placed form, acting on A with a boolean operation
    transformation::Transformation3D{T}          # placement of "right" with respect of "left" 
    left_aabb::AABB{T}
    right_aabb::AABB{T}
end

const BooleanUnion{T} = BooleanShape{T, :union}
const BooleanSubtraction{T} =  BooleanShape{T, :subtraction}
const BooleanIntersection{T} = BooleanShape{T, :intersection}

# Since all the other constructors only use 3 arguments we can sneak the AABB calc in here
function BooleanShape{T,OP}(left,right,place) where {T, OP}
    left_aabb=AABB(extent(left)...)
    right_aabb=AABB(transform_extent(extent(right),place)...)
    BooleanShape{T,OP}(left,right,place,left_aabb,right_aabb)
end

#---Constructor------------------------------------------------------------------------------------
function BooleanUnion(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanUnion{T}(left, right, place)
end
function BooleanSubtraction(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanSubtraction{T}(left, right, place)
end
function BooleanIntersection(left::AbstractShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    BooleanIntersection{T}(left, right, place)
end

function BooleanSubtraction(left::NoShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    return left
end

function BooleanIntersection(left::NoShape{T}, right::AbstractShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    return left
end

function BooleanIntersection(left::AbstractShape{T}, right::NoShape{T}, place::Transformation3D{T}=one(Transformation3D{T})) where T<:AbstractFloat
    return right
end

#---Utilities---------------------------------------------------------------------------------------

#---Printing and Plotting---------------------------------------------------------------------------
function Base.show(io::IO, shape::BooleanUnion{T}) where {T}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanUnion",(left=left, right=right, placement=placement))
end
function Base.show(io::IO, shape::BooleanSubtraction{T}) where {T}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanSubtraction",(left=left, right=right, placement=placement))
end
function Base.show(io::IO, shape::BooleanIntersection{T}) where {T}
    (; left, right, transformation) = shape
    placement = isone(transformation) ? nothing : isone(transformation.rotation) ? transformation.translation : transformation
    print(io, "BooleanIntersection",(left=left, right=right, placement=placement))
end


function GeometryBasics.mesh(shape::BooleanUnion{T}; facets::Integer=64) where {T}
    (; left, right, transformation) = shape
    rmesh = mesh(right; facets)
    merge([mesh(left; facets), Mesh(map(c -> _gbpoint(Point3{T}(c) * transformation), coordinates(rmesh)), faces(rmesh))])
end
function GeometryBasics.mesh(shape::BooleanSubtraction{T}; facets::Integer=64) where {T}
    (; left, right, transformation) = shape
    rmesh = mesh(right; facets)
    merge([mesh(left; facets), Mesh(map(c -> _gbpoint(Point3{T}(c) * transformation), coordinates(rmesh)), faces(rmesh))])
end
function GeometryBasics.mesh(shape::BooleanIntersection{T}; facets::Integer=64) where {T}
    (; left, right, transformation) = shape
    rmesh = mesh(right; facets)
    merge([mesh(left; facets), Mesh(map(c -> _gbpoint(Point3{T}(c) * transformation), coordinates(rmesh)), faces(rmesh))])
end

#---Basic functions---------------------------------------------------------------------------------
function extent(shape::BooleanUnion{T})::Tuple{Point3{T},Point3{T}} where {T}
    (; left, right, transformation) = shape
    minLeft, maxLeft = extent(left)
    minRight, maxRight = transform_extent(extent(right),transformation)
    (min.(minLeft, minRight), max.(maxLeft, maxRight))
end
function extent(shape::BooleanSubtraction{T})::Tuple{Point3{T},Point3{T}} where {T}
    extent(shape.left)
end
function extent(shape::BooleanIntersection{T})::Tuple{Point3{T},Point3{T}} where {T}
     (; left, right, transformation) = shape
    minLeft, maxLeft = extent(left)
    minRight, maxRight = transform_extent(extent(right),transformation)
    (max.(minLeft, minRight), min.(maxLeft, maxRight))
end

 
function inside(shape::BooleanUnion{T}, point::Point3{T})::Int64 where {T}

    course_inside_left = inside(shape.left_aabb, point)
    course_inside_right = inside(shape.right_aabb, point)

    # If the point isn't in either of the BB it is not inside
    !course_inside_left && !course_inside_right && return kOutside
    
    (; left, right, transformation) = shape
    
    # If inside one BB but not the other
    course_inside_left && !course_inside_right && return inside(left, point)
    lpoint = transformation * point
    course_inside_right && !course_inside_left && return inside(right, lpoint)

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

function inside(shape::BooleanIntersection{T}, point::Point3{T})::Int64  where {T}
    
    course_inside_left = inside(shape.left_aabb, point)
    course_inside_right = inside(shape.right_aabb, point)

    # If the point isn't in both BB it is not in the shape
    !(course_inside_left && course_inside_right) && return kOutside

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

function inside(shape::BooleanSubtraction{T}, point::Point3{T})::Int64  where {T}
    course_inside_left = inside(shape.left_aabb,point)
    course_inside_right = inside(shape.right_aabb,point)
    (; left, right, transformation) = shape
    # If the point is the left BB but not the right 
    course_inside_left && !course_inside_right && return inside(left, point)

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

#---Normals------------------------------------------------------------------------------------------
_unit_or_zero(n::Vector3{T}) where T = n ⋅ n > 0 ? normalize(n) : zero(Vector3{T})

function normal(shape::BooleanUnion{T}, point::Point3{T}) where {T}
    (; left, right, transformation) = shape
    lpoint = transformation * point
    posA, posB = inside(left, point), inside(right, lpoint)
    if posA == kSurface && posB == kOutside
        return normal(left, point)
    elseif posB == kSurface && posA == kOutside
        return normal(right, lpoint) * transformation
    elseif posA == kSurface && posB == kSurface      # overlapping faces (touching faces give zero)
        return _unit_or_zero(normal(left, point) + normal(right, lpoint) * transformation)
    end
    return zero(Vector3{T})
end

function normal(shape::BooleanSubtraction{T}, point::Point3{T}) where {T}
    (; left, right, transformation) = shape
    lpoint = transformation * point
    posA, posB = inside(left, point), inside(right, lpoint)
    if posA == kSurface && posB == kOutside
        return normal(left, point)
    elseif posB == kSurface && posA == kInside
        return -(normal(right, lpoint) * transformation)
    elseif posA == kSurface && posB == kSurface      # edge where the cut meets the outer surface
        return _unit_or_zero(normal(left, point) - normal(right, lpoint) * transformation)
    end
    return zero(Vector3{T})
end

function normal(shape::BooleanIntersection{T}, point::Point3{T}) where {T}
    (; left, right, transformation) = shape
    lpoint = transformation * point
    posA, posB = inside(left, point), inside(right, lpoint)
    if posA == kSurface && posB == kInside
        return normal(left, point)
    elseif posB == kSurface && posA == kInside
        return normal(right, lpoint) * transformation
    elseif posA == kSurface && posB == kSurface
        return _unit_or_zero(normal(left, point) + normal(right, lpoint) * transformation)
    end
    return zero(Vector3{T})
end

#---Capacity and surface (no closed form: estimated by sampling, as in Geant4)------------------------
capacity(shape::BooleanShape; npoints::Integer=1_000_000) = estimate_capacity(shape; npoints)
surface(shape::BooleanShape; npoints::Integer=1_000_000) = estimate_surface(shape; npoints)

function safety_to_out(shape::BooleanUnion{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_out(shape::BooleanSubtraction{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_out(shape::BooleanIntersection{T}, point::Point3{T})::T where {T}
    return 0
end

function safety_to_in(shape::BooleanUnion{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_in(shape::BooleanSubtraction{T}, point::Point3{T})::T where {T}
    return 0
end
function safety_to_in(shape::BooleanIntersection{T}, point::Point3{T})::T where {T}
    return 0
end

function distance_to_out(shape::BooleanUnion{T}, point::Point3{T}, dir::Vector3{T})::T where {T}

    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    (; left, right, transformation) = shape
    # If it intersects one BB but not the other just use that one 
    course_intersects_left && !course_intersects_right && return distance_to_out(left, point, dir)
    course_intersects_right && !course_intersects_left && return distance_to_out(right, transformation * point, transformation * dir)

    dist = T(0)
    positionA = inside(left, point)
    if positionA != kOutside  # point inside A
        while(true)
            distA = distance_to_out(left, point, dir)
            dist += (distA > 0 && distA < T(Inf)) ? distA : 0
            dist += push_tolerance(T)
            npoint = point + dist * dir
            lpoint = transformation * npoint
            if inside(right, lpoint) != kOutside # B could be overlapping with A -- and/or connecting A to another part of A
                ldir = transformation * dir
                distB = distance_to_out(right, lpoint, ldir)
                dist += (distB > 0 && distB < T(Inf)) ? distB : 0
                dist += push_tolerance(T)
                npoint = point + dist * dir
                if inside(left, npoint) == kOutside
                    break
                end
            else
                break
            end
        end
        return dist - push_tolerance(T)
    else
        lpoint = transformation * point
        positionB = inside(right, lpoint)
        if positionB != kOutside  # point inside B
            ldir = transformation * dir
            while(true)
                distB = distance_to_out(right, lpoint, ldir)
                dist += (distB > 0 && distB < T(Inf)) ? distB : 0
                dist += push_tolerance(T) # Give a push
                npoint = point + dist * dir
                if inside(left, npoint) != kOutside # A could be overlapping with B -- and/or connecting B to another part of B
                    ldir = transformation * dir
                    distA = distance_to_out(left, npoint, dir)
                    dist += (distA > 0 && distA < T(Inf)) ? distA : 0
                    dist += push_tolerance(T)
                    npoint = point + dist * dir
                    lpoint = transformation * npoint
                    if inside(right, lpoint) == kOutside
                        break
                    end
                else
                    break
                end
            end
            return dist - push_tolerance(T)
        else
            return T(-1)
        end
    end 
end

function distance_to_out(shape::BooleanIntersection{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    (; left, right, transformation) = shape
    distA = distance_to_out(left, point, dir)
    distB = distance_to_out(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distance_to_out(shape::BooleanSubtraction{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    (; left, right, transformation) = shape
    distA = distance_to_out(left, point, dir)
    distB = distance_to_in(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distance_to_in(shape::BooleanUnion{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    # If it doesn't intersect either BB skip everything
    !course_intersects_left && !course_intersects_right && return T(Inf)
    
    (; left, right, transformation) = shape
    
    # If it intersects one BB but not the other just use that one 
    course_intersects_left && !course_intersects_right && return distance_to_in(left, point, dir)
    course_intersects_right && !course_intersects_left && return distance_to_in(right, transformation * point, transformation * dir)
    
    # If it gets here we have to calculate the distance the complicated way
    distA = distance_to_in(left, point, dir)
    distB = distance_to_in(right, transformation * point, transformation * dir)
    return min(distA, distB)
end

function distance_to_in(shape::BooleanIntersection{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    # If it doesn't intersect both BB skip everything
    !(course_intersects_left && course_intersects_right) && return T(Inf)

    # If it gets here we have to calculate the distance the complicated way
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
            d1 = distance_to_in(left, npoint, dir)
            d1 = max(d1, tolerance(T))
            d1 == T(Inf) && return T(Inf)
        end
        if !inRight
            d2 = distance_to_in(right, lpoint, ldir)
            d2 = max(d2, tolerance(T))
            d2 == T(Inf) && return T(Inf)
        end
        if d1 > d2
            # propagate to left shape
            dist += d1
            inleft = true
            npoint += d1 * dir
            lpoint = transformation * npoint
            # check if propagated point is inside right shape, check is done with a little push
            inRight = inside(right, lpoint + tolerance(T) * ldir) == kInside
            inRight && return dist
            # here inleft=true, inright=false
        else
            # propagate to right shape
            dist += d2
            inright = true
            # check if propagated point is inside left shape, check is done with a little push
            npoint += d2 * dir
            lpoint = transformation * npoint
            inLeft = inside(left, npoint + tolerance(T) * dir) == kInside
            inLeft && return dist
        end
    end
    return dist
end

function distance_to_in(shape::BooleanSubtraction{T}, point::Point3{T}, dir::Vector3{T})::T where {T}
    

    invdir=inv.(dir)
    course_intersects_left = intersect(shape.left_aabb, point, dir, invdir)
    course_intersects_right = intersect(shape.right_aabb, point, dir, invdir)

    # If it doesn't intersect the left shape BB we can skip everyingthing
    !course_intersects_left && return T(Inf)

    (; left, right, transformation) = shape

    # If it doesn't intersect the right we can us just the left distance function 
    !course_intersects_right && return distance_to_in(left, point, dir)

    # If it gets here we have to calculate the distance the complicated way
    lpoint = transformation * point
    ldir = transformation * dir
    positionB = inside(left, lpoint)
    inRight = positionB == kInside

    npoint = point
    offset = T(0)   # distance from point to the current probe point npoint
    dist = T(0)

    while true
        if inRight
            # propagate to outside of '- / RightShape'
            d1 = distance_to_out(right, lpoint, ldir)
            dist += (d1 >= 0 && d1 < T(Inf)) ? d1 + push_tolerance(T) : 0
            offset = dist + push_tolerance(T)
            npoint = point + offset * dir
            lpoint = transformation * npoint
            # now master outside 'B'; check if inside 'A'
            inside(left, npoint) == kInside && distance_to_out(left, npoint, dir) > push_tolerance(T) && return dist
        end

        # if outside of both we do a max operation master outside '-' and outside '+' ;  find distances to both
        d2 = distance_to_in(left, npoint, dir)
        d2 = max(d2, 0)
        d2 == T(Inf) && return T(Inf)

        d1 = distance_to_in(right, lpoint, ldir)
        if d2 < d1 - tolerance(T)
          return offset + d2   # d2 is measured from the probe point, not from the original point
        end

        #   propagate to '-'
        dist += (d1 >= 0 && d1 < T(Inf)) ? d1 : 0
        offset = dist + push_tolerance(T)
        npoint = point + offset * dir
        lpoint = transformation * npoint
        inRight = true
    end
end

