#---Material-----------------------------------------------------------------------
struct Isotope{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    N::Int32
    Z::Int32
    Amass::T
end

struct Element{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    composition::Vector{@NamedTuple{fraction::T, isotope::Isotope{T}}}
end

struct Material{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    state::String
    density::T
    temperature::T
    Amass::T
    composition::Vector{@NamedTuple{fraction::T, element::Element{T}}}
    function Material{T}(name::String; state::String="solid", density=0,
                         temperature=293.15, mass=0, composition=[]) where T
        new(name, state, density, temperature, mass, composition)
    end
end



