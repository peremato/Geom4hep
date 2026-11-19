#---Material-----------------------------------------------------------------------
"""
    Isotope{T}(name, N, Z, Amass)

Isotope with `N` nucleons, atomic number `Z` and molar mass `Amass`.
"""
struct Isotope{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    N::Int32
    Z::Int32
    Amass::T
end

"""
    Element{T}(name, composition)

Chemical element: a vector of named tuples `(fraction=..., isotope=...)` with the abundance of each
[`Isotope`](@ref).
"""
struct Element{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    composition::Vector{@NamedTuple{fraction::T, isotope::Isotope{T}}}
end

"""
    Material{T}(name; state="solid", density=0, temperature=293.15, mass=0, composition=[])

Material of a volume. `density` is used by the navigation examples (e.g. [`xray`](@ref)); the GDML
reader stores it in the units of the file (usually g/cm³). `composition` is a vector of named tuples
`(fraction=..., element=...)` with the mass fraction of each [`Element`](@ref).
"""
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



