#---Box-------------------------------------------------------------------------
struct Box{T<:AbstractFloat} <: AbstractShape{T}
    fDimentions::SVector{3,T} # the HALF lengths of the box
    #dx::T
    #dy::T
    #dz::T
end
Box(x,y,z) = Box(SVector(x,y,z))
Base.getindex(b::Box, s::Symbol) = b.fDimentions[coordmap[s]]
Base.getindex(b::Box, i::Int64) = b.fDimentions[i]
#Base.getproperty(b::Box, s::Symbol) = (s in keys(coordmap) ? getfield(b, :fDimentions)[coordmap[s]] : getfield(b,s))
capacity(b::Box) = 8.0 * b[1] * b[2] * b[3]
surface(b::Box) = 8.0(b[1] * b[2] + b[2] * b[3] + b[1] * b[3])
extent(b::Box) = (-b.fDimentions, b.fDimentions)
function normal(box::Box{T}, point::Point3{T}) where T
    safety = abs.(abs.(point) - box.fDimentions)
    if minimum(safety) ≈ 0.0
        n = Vec3{T}(map( t -> t[1] ≈ 0.0 ? sign(t[2]) : 0.0, zip(safety, point)))
        normalize(n)
    else
        nothing    
    end
end
