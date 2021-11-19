#---Material-----------------------------------------------------------------------
struct Isotope <: AbstractMaterial
    name::String
    N::Int32
    Z::Int32
    Amass::Float64
end

struct Element <: AbstractMaterial
    name::String
    composition::Vector{@NamedTuple{fraction::Float64, isotope::Isotope}}
end

struct Material <: AbstractMaterial
    name::String
    state::String
    density::Float64
    temperature::Float64
    Amass::Float64
    composition::Vector{@NamedTuple{fraction::Float64, element::Element}}
    function Material(name::String; state::String="solid", density::Float64=0.0,
                      temperature::Float64=293.15, mass::Float64=0.0, composition=[])
        new(name, state, density, temperature, mass, composition)
    end
end



