using StaticArrays
using BenchmarkTools
using Virtual

const Point3 = SVector{3}
const Vector3 = SVector{3}

abstract type AbstractShape{T<:AbstractFloat} end

struct Box{T<:AbstractFloat} <: AbstractShape{T}
    fDimensions::Vector3{T}
    Box{T}(x,y,z) where T = new((x,y,z))
end
struct Tube{T<:AbstractFloat} <: AbstractShape{T}
    rmin::T
    rmax::T
    z::T
    Tube{T}(rmin, rmax, z) where T = new(rmin, rmax, z)
end
struct Cone{T<:AbstractFloat} <: AbstractShape{T}
    r::T
    z::T
    Cone{T}(r, z) where T = new(r, z)
end
struct BooleanUnion{T<:AbstractFloat, SL<:AbstractShape{T}, SR<:AbstractShape{T}} <: AbstractShape{T}
    left::SL 
    right::SR
end
BooleanUnion(left::AbstractShape{T}, right::AbstractShape{T}) where T = BooleanUnion{T,typeof(left),typeof(right)}(left, right)


@virtual fast_inside(s::AbstractShape{Float64}, p::Point3{Float64}) = error("No default method for inside!")
@virtual fast_inside(s::AbstractShape{Float32}, p::Point3{Float32}) = error("No default method for inside!")

@override fast_inside(box::Box{Float64}, point::Point3{Float64}) = inside(box, point)
@override fast_inside(box::Box{Float32}, point::Point3{Float32}) = inside(box, point)
@override fast_inside(tube::Tube{Float64}, point::Point3{Float64}) = inside(tube, point)
@override fast_inside(tube::Tube{Float32}, point::Point3{Float32}) = inside(tube, point)
@override fast_inside(cone::Cone{Float64}, point::Point3{Float64}) = inside(cone, point)
@override fast_inside(cone::Cone{Float32}, point::Point3{Float32}) = inside(cone, point)
@override fast_inside(union::BooleanUnion{Float64, SL, SR} where {SL,SR}, point::Point3{Float64}) = inside(union, point)
@override fast_inside(union::BooleanUnion{Float32, SL, SR} where {SL,SR}, point::Point3{Float32}) = inside(union, point)

#---The real algorithms with papametric methods (ideally I need to keep as they are)
inside(box::Box{T}, point::Point3{T}) where T = maximum(abs.(point) - box.fDimensions) < T(0)
inside(tube::Tube{T}, point::Point3{T}) where T = tube.rmin^2 < point[1]^2+point[2]^2 < tube.rmax^2  && abs(point[3]) < tube.z
inside(cone::Cone{T}, point::Point3{T}) where T = point[1]^2+point[2]^2 <  ((cone.z-point[3])/2cone.z)^2 && abs(point[3]) < cone.z
inside(union::BooleanUnion{T,SL,SR}, point::Point3{T}) where {T,SL,SR} = inside(union.left, point) && inside(union.right, point)

box1 = Box{Float64}(1,2,3)
box2 = Box{Float64}(2,2,2)
tube = Tube{Float64}(1,2,3)
cone = Cone{Float64}(2,3)

const samples = (box1, box2, tube, cone ) #, BooleanUnion(box1,tube), BooleanUnion(box2,tube), BooleanUnion(box1,box2) )
shapes = AbstractShape{Float64}[samples[rand(1:length(samples))] for i = 1:100]

const Shape{T} = Union{Box{T}, Tube{T}, Cone{T} #=, BooleanUnion{T}=# } where T<:AbstractFloat
ushapes = Shape{Float64}[samples[rand(1:length(samples))] for i = 1:100]

function inside_score(inside_func, point::Point3{T}, shapes::Vector{S}) where {T,S}
    s = 0
    for shape in shapes
        inside_func(shape, point) && ( s += 1)
    end
    return s
end

point = Point3{Float64}(0.9,1,1)

@info "fast_inside via Virtual.jl"
@btime inside_score(fast_inside, point, shapes)
@info "original inside with Vector of AbstractShapes"
@btime inside_score(inside, point, shapes)
@info "original inside with Vector of union Shape"
@btime inside_score(inside, point, ushapes)

