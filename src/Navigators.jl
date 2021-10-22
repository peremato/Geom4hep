
abstract type AbstractNavigator end
abstract type AbstractNavigatorState end

struct SimpleNavigator{T<:AbstractFloat} <: AbstractNavigator
    top::Volume
end

mutable struct NavigatorState{T<:AbstractFloat} <: AbstractNavigatorState
    top::Volume{T}                               # Typically the unplaced world
    currvol::Volume{T}                           # Current volume
    #volstack::Vector{AbstractPlacedVolume{T}}    # Stack of daughters 
    #tolocal::Transformation3D{T}                 # Composed transform of daughters
    volstack::Vector{Int64}                      # keep the indexes of all daughters up to the current one 
    tolocal::Vector{Transformation3D{T}}         # Stack of transformations
    function NavigatorState{T}(top::Volume{T}) where T<:AbstractFloat
        x = new{T}(top, top, Vector{Int64}(), Vector{Transformation3D{T}}()) 
        sizehint!(x.volstack,16)
        sizehint!(x.tolocal,16)
        x
    end
end

#NavigatorState(top::Volume{T}) where T<:AbstractFloat =  NavigatorState{T}(top, Vector{AbstractPlacedVolume{T}}(), Vector{Transformation3D{T}}())
currentVolume(state::NavigatorState) = state.currvol

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
    step::T = step_limit
    candidate::Int64 = 0
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

function getDaughter(vol::Volume{T}, path::Vector{Int64}) where T<:AbstractFloat
    for idx in path
        vol = vol.daughters[idx].volume
    end
    return vol
end

#---Update the NavigatorState and get the step (distance)
function computeStep!(state::NavigatorState{T}, gpoint::Point3{T}, gdir::Vector3{T}, step_limit::T) where T<:AbstractFloat
    #lpoint, ldir = isone(state.tolocal) ? (gpoint, gdir) : (state.tolocal * gpoint, state.tolocal * gdir)
    lpoint = transform(state.tolocal, gpoint)
    ldir = transform(state.tolocal, gdir)
    volume = state.currvol

    step, idx = getClosestDaughter(volume, lpoint, ldir, step_limit)

    #---If didn't hit any daughter return distance to out
    if idx == 0
        step = distanceToOut(volume.shape, lpoint, ldir)
        step = step < 0. ? 0. : step
        if step > 0.
            if !isempty(state.volstack) 
                pop!(state.volstack)
                pop!(state.tolocal)
                state.currvol = getDaughter(state.top, state.volstack)
            end
        end
    end
    #---We hit a daughter, push it into the stack
    if idx != 0
        step += kTolerance # to ensure that do not stay in the surface of the daughter
        pvol = volume.daughters[idx]
        push!(state.volstack, idx)
        push!(state.tolocal, pvol.transformation)
        state.currvol = pvol.volume   
    end
    return step
end
