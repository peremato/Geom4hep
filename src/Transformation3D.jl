#---Transformation3D----------------------------------------------------------------------
struct Transformation3D{T<:AbstractFloat}
    rotation::RotMatrix3{T}
    translation::Vector3{T}
end

Transformation3D() = Transformation3D{Float64}()

# Usefull constructors 
Transformation3D{T}(dx, dy, dz) where T<:AbstractFloat = Transformation3D(one(RotMatrix3{T}), Vector3{T}(dx, dy, dz))
Transformation3D{T}(dx, dy, dz, rotx, roty, rotz) where T<:AbstractFloat = Transformation3D(RotMatrix3{T}(RotXYZ{T}(rotx, roty, rotz)), Vector3{T}(dx, dy, dz))
Transformation3D{T}(dx, dy, dz, rot::Rotation{3,T}) where T<:AbstractFloat = Transformation3D(RotMatrix3{T}(rot), Vector3{T}(dx, dy, dz))

# Transforms
@inline transform(t::Transformation3D{T}, p::Point3{T}) where T<:AbstractFloat = t.rotation * (p - t.translation)
@inline transform(t::Transformation3D{T}, d::Vector3{T}) where T<:AbstractFloat = t.rotation * d
@inline invtransform(t::Transformation3D{T}, p::Point3{T}) where T<:AbstractFloat = t.translation + (p' * t.rotation)'
@inline invtransform(t::Transformation3D{T}, d::Vector3{T}) where T<:AbstractFloat = (d' * t.rotation)'
@inline Base.:*(t::Transformation3D{T}, p::Point3{T}) where T<:AbstractFloat =  transform(t,p)
@inline Base.:*(t::Transformation3D{T}, d::Vector3{T}) where T<:AbstractFloat =  transform(t,d)
@inline Base.:*(p::Point3{T}, t::Transformation3D{T}) where T<:AbstractFloat =  invtransform(t,p)
@inline Base.:*(d::Vector3{T}, t::Transformation3D{T}) where T<:AbstractFloat =  invtransform(t,d)

function transform(v::Vector{Transformation3D{T}}, p::Point3{T}) where T<:AbstractFloat
    for t in v
        p = t.rotation * (p - t.translation)
    end
    return p
end
function transform(v::Vector{Transformation3D{T}}, d::Vector3{T}) where T<:AbstractFloat
    for t in v
        d = t.rotation * d
    end
    return d
end
function transform(v::Vector{Transformation3D{T}}, p::Point3{T}, d::Vector3{T}) where T<:AbstractFloat
    for t in v
        p = t.rotation * (p - t.translation)
        d = t.rotation * d
    end
    return p, d
end

# Compose
Base.:*(t1::Transformation3D{T}, t2::Transformation3D{T}) where T<:AbstractFloat = 
    Transformation3D{T}(t1.rotation * t2.rotation, t1.translation + t1.rotation * t2.translation)

# Utilities
Base.one(::Type{Transformation3D{T}}) where T<:AbstractFloat = Transformation3D{T}(one(RotMatrix3{T}), Vector3{T}(0,0,0))
Base.isone(t::Transformation3D{T}) where T<:AbstractFloat = isone(t.rotation) && iszero(t.translation)
hasrotation(t::Transformation3D{T}) where T<:AbstractFloat = !isone(t.rotation)
hastranslation(t::Transformation3D{T}) where T<:AbstractFloat = !iszero(t.translation)

