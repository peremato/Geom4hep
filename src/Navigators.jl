
abstract type AbstractNavigator end
abstract type AbstractNavigatorState end

struct SimpleNavigator{T<:AbstractFloat} <: AbstractNavigator
    top::Volume{T}
end

mutable struct NavigatorState{T<:AbstractFloat} <: AbstractNavigatorState
    top::Volume{T}                               # Typically the unplaced world
    volstack::Vector{Int64}                      # keep the indexes of all daughters up to the current one 
    tolocal::Vector{Transformation3D{T}}         # Stack of transformations
    function NavigatorState{T}(top::Volume{T}) where T<:AbstractFloat
        x = new{T}(top, Vector{Int64}(), Vector{Transformation3D{T}}()) 
        sizehint!(x.volstack,16)
        sizehint!(x.tolocal,16)
        x
    end
end

function reset!(state::NavigatorState{T}) where T<:AbstractFloat
    empty!(state.volstack)
    empty!(state.tolocal)
end

function currentVolume(state::NavigatorState)
    vol = state.top
    for idx in state.volstack
        vol = vol.daughters[idx].volume
    end
    return vol
end

function collectDaughters!(path::Vector{Int64}, transforms::Vector{Transformation3D{T}}, vol::Volume{T}, point::Point3{T}) where T<:AbstractFloat
    for (idx, daughter) in enumerate(vol.daughters)
        if contains(daughter, point)
            push!(path, idx)
            push!(transforms, daughter.transformation)
            collectDaughters!(path, transforms, daughter.volume, daughter.transformation * point)
            break
        end
    end
end

#---Locate global point and initialize state 
function locateGlobalPoint!(state::NavigatorState{T}, gpoint::Point3{T}) where T<:AbstractFloat
    vol = state.top
    empty!(state.volstack)
    empty!(state.tolocal)
    if contains(vol, gpoint)
        collectDaughters!(state.volstack, state.tolocal, vol, gpoint)
        return true
    end
    return false
end

function getClosestDaughter(volume::Volume{T}, point::Point3{T}, dir::Vector3{T}, step_limit::T ) where T<:AbstractFloat
    step = step_limit
    candidate = 0
    #---Linear loop over the daughters
    for (idx, daughter) in enumerate(volume.daughters)
        if !contains(daughter, point)
            dist = distanceToIn(daughter, point, dir)
            if dist < step && dist != Inf && dist > 0.
                step = dist
                candidate = idx
            end
        else
            step = -1.
        end
    end
    return step, candidate
end

#---Update the NavigatorState and get the step (distance)
function computeStep!(state::NavigatorState{T}, gpoint::Point3{T}, gdir::Vector3{T}, step_limit::T) where T<:AbstractFloat
    lpoint, ldir = transform(state.tolocal, gpoint, gdir)
    #volume = state.currvol
    volume = currentVolume(state)

    step, idx = getClosestDaughter(volume, lpoint, ldir, step_limit)

    #---If didn't hit any daughter return distance to out
    if idx == 0
        step = distanceToOut(volume.shape, lpoint, ldir)
        step = step < 0. ? 0. : step
        if step > 0.
            if !isempty(state.volstack) 
                pop!(state.volstack)
                pop!(state.tolocal)
            end
        end
    end
    #---We hit a daughter, push it into the stack
    if idx != 0
        step += kTolerance # to ensure that do not stay in the surface of the daughter
        pvol = volume.daughters[idx]
        push!(state.volstack, idx)
        push!(state.tolocal, pvol.transformation)
    end
    return step
end
