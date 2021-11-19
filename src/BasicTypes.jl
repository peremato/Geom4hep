#--- basic stuff

abstract type AbstractShape{T<:AbstractFloat} end
abstract type AbstractMaterial end

const Vector3 = SVector{3}
Base.:*(x::Vector3, y::Vector3) = Vector3(x[1]*y[1], x[2]*y[2], x[3]*y[3])
Base.:*(x::Vector3, n::Number) = Vector3(x[1]*n, x[2]*n, x[3]*n)
Base.:*(n::Number, x::Vector3) = x * n
Base.:/(p::Point3, n::Number) = Point3(p[1]/n, p[2]/n, p[3]/n)
LinearAlgebra.:⋅(x::Vector3, y::Vector3) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
Base.:-(p1::Point3, p2::Point3) = Vector3(p1[1]-p2[1], p1[2]-p2[2], p1[3]-p2[3])
Base.:+(p1::Point3, p2::Point3) = Vector3(p1[1]+p2[1], p1[2]+p2[2], p1[3]+p2[3])
LinearAlgebra.:×(x::Vector3, y::Vector3) = Vector3(x[2]*y[3]-x[3]*y[2], -x[1]*y[3]+x[3]*y[1], x[1]*y[2]-x[2]*y[1])
unitize(v::Vector3) = v/(v⋅v)

#---constants
const kTolerance = 1e-9
const kInside  = 0
const kSurface = 1
const kOutside = 2
export kInside, kSurface, kOutside

#=
#---dymanic (mutable) points and vectors
mutable struct DPoint3{T} <: AbstractVector{T}
    x::T
    y::T
    z::T
end

function move!(p::DPoint3{T}, d::Vector3{T}, s::T) where T<:AbstractFloat
    p .=  p .+ d .* s
end

#---Array interface for DPoint
Base.size(::DPoint3{T}) where T = (3,)
function Base.getindex(p::DPoint3{T}, i::Int) where T
    if i == 1;     return p.x
    elseif i == 2; return p.y
    elseif i == 3; return p.z
    else    
        throw(BoundsError(p, i))
    end   
end
function Base.setindex!(p::DPoint3{T}, v::N, i::Int) where {T,N}
    if i == 1;     p.x = v
    elseif i == 2; p.y = v
    elseif i == 3; p.z = v
    else    
      throw(BoundsError(p, i))
    end 
end
Base.convert(::Type{Point3{T}}, p::DPoint3{T}) where T = Point3{T}(p)
=#





