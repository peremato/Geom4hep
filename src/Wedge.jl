#---Wedge -------------------------------------------------------------
struct Wedge{T<:AbstractFloat} <: AbstractShape{T}
    ϕ₀::T             # starting angle
    Δϕ::T             # delta angle representing/defining the wedge
    along1::Vector2{T}  # vector along the first plane
    along2::Vector2{T}  # vector along the second plane
    normal1::Vector2{T} # normal vector for first plane
    normal2::Vector2{T} # normal vector for second plane
    function Wedge{T}(ϕ₀, Δϕ) where T<:AbstractFloat
        along1 = Vector2{T}(cos(ϕ₀), sin(ϕ₀))
        along2 = Vector2{T}(cos(ϕ₀+Δϕ), sin(ϕ₀+Δϕ))
        normal1 = Vector2{T}(-sin(ϕ₀), cos(ϕ₀))
        normal2 = Vector2{T}(sin(ϕ₀+Δϕ), -cos(ϕ₀+Δϕ))
        new(ϕ₀, Δϕ, along1, along2, normal1, normal2)
    end
end

function Base.show(io::IO, wed::Wedge{T}) where T
    print(io, "Wedge{$T}",(wed.ϕ₀,wed.Δϕ))
end

function isOnSurface(along::Vector2{T}, normal::Vector2{T}, x::T, y::T) where T<:AbstractFloat  
    along[1]*x + along[2]*y >= 0. && abs(normal[1]*x + normal[2]*y) < kTolerance(T)
end

function isInside(w::Wedge{T}, x::T, y::T; tol::T=T(0)) where T<:AbstractFloat
    w.Δϕ == T(2π) && return true
    startCheck = (-x * w.along1[2] + y * w.along1[1])
    endCheck   = (-w.along2[1] * y + w.along2[2] * x)
    if w.Δϕ < T(π)
      inside = startCheck > -tol && endCheck > -tol
    else
      inside = startCheck > -tol || endCheck > -tol
    end
    inside &= !isOnSurface(w.along1, w.normal1, x, y) && !isOnSurface(w.along2, w.normal2, x, y)
end

function isOutside(w::Wedge{T}, x::T, y::T; tol::T=T(0) ) where T<:AbstractFloat
    w.Δϕ == T(2π) && return false
    startCheck = (-x * w.along1[2] + y * w.along1[1])
    endCheck   = (-w.along2[1] * y + w.along2[2] * x)
    if w.Δϕ < T(π)
      outside = startCheck < tol || endCheck < tol
    else
      outside = startCheck < tol && endCheck < tol
    end
    outside &= !isOnSurface(w.along1, w.normal1, x, y) && !isOnSurface(w.along2, w.normal2, x, y)
end
isOutside(w::Wedge{T}, p::Point2{T}; tol::T=T(0) ) where T<:AbstractFloat = isOutside(w, p[1], p[2]; tol=tol)
isInside(w::Wedge{T}, p::Point2{T}; tol::T=T(0) ) where T<:AbstractFloat = isInside(w, p[1], p[2]; tol=tol)

@override function safetyToOut(w::Wedge{T}, x::T, y::T) where T<:AbstractFloat
    (; normal1, normal2, Δϕ) = w
    dist1 = x * normal1[1] + y * normal1[2]
    dist2 = x * normal2[1] + y * normal2[2]
    if Δϕ < π 
      min(dist1, dist2)
    else
      max(dist1, dist2)
    end
end