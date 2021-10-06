
abstract type AbstractNavigator end
abstract type AbstractNavigatorState end

struct SimpleNavigator{T<:AbstractFloat} <: AbstractNavigator
    top::Volume
end

mutable struct NavigatorState{T<:AbstractFloat} <: AbstractNavigatorState
    top::Volume{T}                               # Typically the unplaced world
    volstack::Vector{AbstractPlacedVolume{T}}    # Stack of daughters 
    tolocal::Transformation3D{T}                 # Composed transform of daughters
end

NavigatorState(top::Volume{T}) where T<:AbstractFloat =  NavigatorState{T}(top, Vector{AbstractPlacedVolume{T}}(), one(Transformation3D{T}))
currentVolume(state::NavigatorState) = isempty(state.volstack) ? state.top : last(state.volstack).volume

function collectDaughters!(path::Vector{AbstractPlacedVolume{T}}, vol::Volume{T}, point::Point3{T}) where T<:AbstractFloat
    for daughter in vol.daughters
        if contains(daughter, point)
            push!(path, daughter)
            collectDaughters!(path, daughter.volume, daughter.transformation * point)
            break
        end
    end
end

#---Locate global point and initialize state 
function locateGlobalPoint!(state::NavigatorState, vol::Volume{T}, gpoint::Point3{T}) where T<:AbstractFloat
    state.top = vol
    state.volstack = Vector{AbstractPlacedVolume{T}}()
    state.tolocal = one(Transformation3D{T})
    if contains(vol, gpoint)
        collectDaughters!(state.volstack, vol, gpoint)
        state.tolocal = prod(x -> getfield(x,:transformation), state.volstack, init = one(Transformation3D{T}))
        return true
    end
    return false
end

#---Update the NavigatorState and get the step (distance)
function computeStep!(state::NavigatorState, gpoint::Point3{T}, gdir::Vector3{T}, step_limit::T) where T<:AbstractFloat
    lpoint, ldir = isone(state.tolocal) ? (gpoint, gdir) : (state.tolocal * gpoint, state.tolocal * gdir)
    volume = isempty(state.volstack) ? state.top : last(state.volstack).volume
    step      = step_limit
    candidate = nothing

    #---Linear loop over the daughters 
    for daughter in volume.daughters
        if !contains(daughter, lpoint)
            dist = distanceToIn(daughter, lpoint, ldir)
            valid = dist < step && dist != Inf && dist > 0.
            step = valid ? dist : step
            candidate = valid ? daughter : candidate
        else
            step = -1.
        end
    end

    #---If didn't hit any daughter return distance to out
    if candidate == nothing
        step = distanceToOut(volume.shape, lpoint, ldir)
        step = step < 0. ? 0. : step
        if step > 0.
            if !isempty(state.volstack) 
                pop!(state.volstack)
            end
            state.tolocal = prod(x -> getfield(x,:transformation), state.volstack, init = one(Transformation3D{T}))
        end
    end

    if candidate != nothing
        step += kTolerance # to ensure that do not stay in the surface of the daughter
        push!(state.volstack, candidate)
        state.tolocal = candidate.transformation * state.tolocal
    end
    return step
end
