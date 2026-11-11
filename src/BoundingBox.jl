#---Bounding Boxes---------------------------------------------------------------------------------
# Axis Aligned Bounding Box, defined by the min/max coordinates
struct AABB{T<:AbstractFloat}
    min::Point3{T}
    max::Point3{T}
end

_area(d) = 2(d[1]*d[2] + d[2]*d[3] + d[3]*d[1])
_volume(d) = d[1]*d[2]*d[3]

diagonal(bb::AABB{T}) where T = bb.max .- bb.min
area(bb::AABB{T}) where T = _area(diagonal(bb))
volume(bb::AABB{T}) where T = _volume(diagonal(bb))

function GeometryBasics.coordinates(bb::AABB{T}) where T
    a, b = bb.min, bb.max
    Point3{T}[a, (b[1], a[2], a[3]), (a[1], b[2], a[3]), (a[1], a[2], b[3]),
                 (a[1], b[2], b[3]), (b[1], a[2], b[3]), (b[1], b[2], a[3]), b]
end
GeometryBasics.faces(::AABB{T}) where T = QuadFace{Int64}[(1,3,7,2), (1,2,6,4), (1,4,5,3), (8,6,2,7), (8,7,3,5), (8,5,4,6)]
GeometryBasics.normals(::AABB{T}) where T = Point3{T}[(0,0,-1), (0,-1,0), (-1,0,0), (1,0,0), (0,1,0), (0,0,1)]

function intersect(bb::AABB{T}, point::Point3{T}, dir::Vector3{T}, rcp_dir::Vector3{T}=inv.(dir)) where {T}
    (distance, distout) =  intersect_aabox_ray(bb.min, bb.max, point, dir, rcp_dir)
    (distance >= distout  || distout <= tolerance(T) / 2 ) ? false : true
end

function inside(b::AABB{T}, point::Point3{T}) where {T}
    all(b.min .<= point) && all(b.max .>= point)
end

#---Fast Axial Aligned Bounded Box Ray intersection routine.  Returns both distances to the intesections
function intersect_aabox_ray(bbmin::Point3{T}, bbmax::Point3{T}, point::Point3{T}, dir::Vector3{T}, rcp_dir::Vector3{T}=inv.(dir)) where {T}
    
    tmin = typemin(T)
    tmax = typemax(T)

    t1v = (bbmin - point) * rcp_dir
    t2v = (bbmax - point) * rcp_dir
    
    # From here on down all we are doing is calculating the following commented lines
    # Importantly this is non branching and fast without using fastmath 
    # tmin = maximum(min.(t1v,t2v))
    # tmax = minimum(max.(t1v,t2v))

    min(x, y) = ifelse(x < y, x, y)
    max(x, y) = ifelse(x > y, x, y)


    for i in 1:3
        t1 = t1v[i]
        t2 = t2v[i]
        flip = t1 > t2
        tmin = max(ifelse(flip, t2, t1), tmin)
        tmax = min(ifelse(flip, t1, t2), tmax)
    end
    return (tmin,tmax)
end

function transform_extent(extent::Tuple{Point3{T},Point3{T}}, transformation) where {T}
    a, b = extent
    corners = (a, Point3{T}(b[1], a[2], a[3]), Point3{T}(a[1], b[2], a[3]), Point3{T}(a[1], a[2], b[3]),
                  Point3{T}(a[1], b[2], b[3]), Point3{T}(b[1], a[2], b[3]), Point3{T}(b[1], b[2], a[3]), b)
    pmin = Point3{T}(typemax(T), typemax(T), typemax(T))
    pmax = Point3{T}(typemin(T), typemin(T), typemin(T))
    for c in corners
        p = c * transformation
        pmin = min.(pmin, p)
        pmax = max.(pmax, p)
    end
    return (pmin, pmax)
end