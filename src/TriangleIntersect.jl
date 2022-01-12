struct Triangle{T<:AbstractFloat}
    a::Point3{T}
    v1::Vector3{T}
    v2::Vector3{T}
    normal::Vector3{T}
    v1v1::T
    v2v2::T
    v1v2::T
    denom::T
    function Triangle{T}(a::Point3{T}, b::Point3{T}, c::Point3{T}) where T <: AbstractFloat
        v1 = b-a
        v2 = c-a
        normal = v1 × v2
        v1v1 = v1⋅v1
        v2v2 = v2⋅v2
        v1v2 = v1⋅v2
        denom = v1v2*v1v2 - v1v1*v2v2
        new(a, v1, v2, normal, v1v1, v2v2, v1v2, denom)
    end
end

struct Intersection
    ip::Point3{Float64} # intersecting point
    id::Float64 # intersecting distance
    is_intersection::Bool
end

const no_intersection = (0.0, false)

function Base.intersect(p::Point3, d::Vector3, t::Triangle; out::Bool=true)
    # denom = t.normal ⋅ d
    n1, n2, n3 = t.normal
    a1, a2, a3 = t.a
    p1, p2, p3 = p
    d1, d2, d3 = d
    denom = d1*n1 + d2*n2 + d3*n3
    denom == 0 && return no_intersection
    #out && denom < 0 && return no_intersection

    #ri = t.normal⋅(t.a - p) / denom
    ri =  (n1*(a1 - p1) + n2*(a2 - p2) + n3*(a3 - p3)) / denom
    ri <= 0 && return no_intersection

    #plane_intersection =  ri * d + p
    #w = plane_intersection - t.a
    #wv1 = w ⋅ t.v1
    #wv2 = w ⋅ t.v2
    w1 = ri * d1 + p1 - a1
    w2 = ri * d2 + p2 - a2
    w3 = ri * d3 + p3 - a3
    wv1 = w1 * t.v1[1] + w2 * t.v1[2] + w3 * t.v1[3]
    wv2 = w1 * t.v2[1] + w2 * t.v2[2] + w3 * t.v2[3]

    s_intersection = (t.v1v2 * wv2 - t.v2v2 * wv1) / t.denom
    s_intersection <= -kTolerance(T)/2 && return no_intersection
    s_intersection >= 1 + kTolerance(T)/2 && return no_intersection
    t_intersection = (t.v1v2 * wv1 - t.v1v1 * wv2) / t.denom
    t_intersection <= -kTolerance(T)/2 && return no_intersection
    t_intersection >= 1 + kTolerance(T)/2 && return no_intersection
    s_intersection + t_intersection >= 1 + kTolerance(T)/2 && return no_intersection
    #Intersection(t.a + s_intersection*t.v1 + t_intersection*t.v2, ri, true)
    return (ri, true)
end

#=
function distanceToPlane(p::PV, d::PV, a::PV, b::PV, c::PV)
    #a1, a2, a3 = a
    #b1, b2, b3 = b
    #c1, c2, c3 = c
    #p1, p2, p3 = point
    #d1, d2, d3 = direction
    n1 = (b.y-a.y)*(c.z-a.z) - (b.z-a.z)*(c.y-a.y)
    n2 = (b.z-a.z)*(c.x-a.x) - (b.x-a.x)*(c.z-a.z)
    n3 = (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x)
    ((a.x - p.x)*n1 + (a.y - p.y)*n2 + (a.z - p.z)*n3) / (d.x*n1 + d.y*n2 + d.z*n3)
    #n = (b-a) × (c-a)
    #d = ((a - point)⋅n)/(direction⋅n)
end
=#