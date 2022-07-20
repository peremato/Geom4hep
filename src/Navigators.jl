
abstract type AbstractNavigator end
abstract type AbstractNavigatorState end

struct SimpleNavigator{T<:AbstractFloat} <: AbstractNavigator
    top::Volume{T}
end

mutable struct NavigatorState{T<:AbstractFloat} <: AbstractNavigatorState
    topvol::Volume{T}                            # Typically the unplaced world
    currvol::Volume{T}                           # the current volume
    isinworld::Bool                              # inside world volume
    volstack::Vector{Int64}                      # keep the indexes of all daughters up to the current one 
    tolocal::Vector{Transformation3D{T}}         # Stack of transformations
    function NavigatorState{T}(top::Volume{T}) where T<:AbstractFloat
        x = new{T}(top, top, false, Vector{Int64}(), Vector{Transformation3D{T}}()) 
        sizehint!(x.volstack,16)
        sizehint!(x.tolocal,16)
        return x
    end
end

@inline function reset!(state::NavigatorState{T}) where T<:AbstractFloat
    empty!(state.volstack)
    empty!(state.tolocal)
    state.currvol = state.topvol
    state.isinworld = false 
end

@inline function currentVolume(state::NavigatorState{T}) where T<:AbstractFloat
    return state.currvol
end

@inline function popOut!(state::NavigatorState{T}) where T<:AbstractFloat
    if !isempty(state.volstack)
        pop!(state.volstack)
        pop!(state.tolocal)
        vol = state.topvol
        for idx in state.volstack
            vol = vol.daughters[idx].volume
        end
        state.currvol = vol
    else
        state.currvol = state.topvol
        state.isinworld = false
    end
end

@inline function pushIn!(state::NavigatorState{T}, pvol::PlacedVolume{T}) where T<:AbstractFloat
    push!(state.volstack, pvol.idx)
    push!(state.tolocal, pvol.transformation)
    state.currvol = pvol.volume
end

#---Locate global point and initialize state 
function locateGlobalPoint!(state::NavigatorState{T}, point::Point3{T}) where T<:AbstractFloat
    reset!(state)
    vol = state.topvol
    if contains(vol, point)
        state.currvol = vol
        state.isinworld = true
    end
    # check the daughters in recursive manner
    isinside = true
    while isinside
        isinside = false
        for pvol in vol.daughters
            if contains(pvol, point)
                pushIn!(state, pvol)
                isinside = true
                vol   = pvol.volume
                point = pvol.transformation * point
                break
            end
        end
    end
    return state.isinworld
end

@inline function isInVolume(state::NavigatorState{T}) where T<:AbstractFloat
    state.isinworld
end

function getClosestDaughter(volume::Volume{T}, point::Point3{T}, dir::Vector3{T}, step_limit::T ) where T<:AbstractFloat
    step = step_limit
    candidate = 0
    #---Linear loop over the daughters
    for (idx, daughter) in enumerate(volume.daughters)
        #---Assuming that it is not yet inside the daughter (otherwise it returns -1.)
        dist = distanceToIn(daughter, point, dir)
        dist < 0. && return (zero(T), idx) 
        if dist > 0. && dist != Inf && dist < step
            step = dist
            candidate = idx
        end
    end
    return step, candidate
end

#---Update the NavigatorState and get the step (distance)
function computeStep!(state::NavigatorState{T}, gpoint::Point3{T}, gdir::Vector3{T}, step_limit::T) where T<:AbstractFloat
    lpoint, ldir = transform(state.tolocal, gpoint, gdir) 
    volume = currentVolume(state)
    step, idx = getClosestDaughter(volume, lpoint, ldir, step_limit)
    #---If didn't hit any daughter return distance to out
    if idx == 0
        step = distanceToOut(volume.shape, lpoint, ldir)
        if step >= 0.
            popOut!(state)
            step += kTolerance(T)
        end
    else
    #---We hit a daughter, push it into the stack
        step += kTolerance(T) # to ensure that do not stay in the surface of the daughter
        pvol = volume.daughters[idx]
        pushIn!(state, pvol)
    end
    return step
end

#=
("====>", step, lpoint, volume.shape) = ("====>", -8.526512829121202e-14, 
[98.61111111111111, 138.05555555555554, 389.55000000000007], 
Polycone{Float64, 11}(rmin = [22.700000000000003, 22.700000000000003, 22.700000000000003, 22.700000000000003, 22.700000000000003, 123.80000000000001, 123.80000000000001, 22.700000000000003, 22.700000000000003, 22.700000000000003, 22.700000000000003, 22.700000000000003], 
                     rmax = [169.0, 169.0, 177.5, 177.5, 177.5, 177.5, 177.5, 177.5, 177.5, 177.5, 169.0, 169.0], 
                     z = [-398.0, -389.55, -389.55, -317.70000000000005, -296.1, -296.1, 296.1, 296.1, 317.70000000000005, 389.55, 389.55, 398.0], 
                     ϕ₀ = 0.0, Δϕ = 6.283185307179586))
=#