#---Used for Aggregates------------------------------------------------------------- 
"""
    NoShape{T}()

Placeholder for a shape that is not supported (e.g. read from GDML), and the shape of an
[`Assembly`](@ref). It contains no points.
"""
struct NoShape{T,PV} <: AbstractShape{T}
    pvolumes::Vector{PV}
end
NoShape{T}() where T = NoShape{T,Nothing}([])

