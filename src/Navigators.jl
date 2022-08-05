
abstract type AbstractNavigator end
abstract type AbstractNavigatorState end

#---Global navigator--------------------------------------------------------------------------------
#gNavigator = nothing

#---Define the various navigators-------------------------------------------------------------------
struct TrivialNavigator{T<:AbstractFloat} <: AbstractNavigator
    world::Volume{T}
    function TrivialNavigator{T}(world::Volume{T}) where T
        #global gNavigator
        nav = new(world)
        #gNavigator = nav
    end    
end

#---BVHNavigator-------------------------------------------------------------------------------------
struct BVHNavigator{T<:AbstractFloat} <: AbstractNavigator
    world::Volume{T}
    param::BVHParam{T}
    bvhdict::Dict{UInt64,BVH{T}}
end
function BVHNavigator{T}(world::Volume{T}) where T
    global gNavigator
    nav = BVHNavigator{T}(world, BVHParam{T}(8), Dict{UInt64,BVH{T}}())
    #--- Create accelerationstructures---------------------------------
    buildBVH(nav, world, Set{UInt64}())
    #gNavigator =  nav
    nav
end
function buildBVH(nav::BVHNavigator{T}, vol::Volume{T}, set::Set{UInt64}) where T
    id = objectid(vol)
    #---return immediatelly if BVH is already there----
    id in set && return
    push!(set, id)
    #---only for more than 16 daughters
    if length(vol.daughters) > 16
        nav.bvhdict[id] = buildBVH(vol.daughters, nav.param)
    end
    #---Call recursively itself------------------------
    for d in vol.daughters
        buildBVH(nav, d.volume, set)
    end
end

function Base.show(io::IO, nav::BVHNavigator{T}) where T
    v = length(nav.bvhdict)
    b = count(i->isnothing(i), values(nav.bvhdict))
    print(io, "BVHNavigator{$T}",(world=nav.world, logicalvols=v, bvhstructs=(v - b)))
end

#---NavigatorState----------------------------------------------------------------------------------
mutable struct NavigatorState{T<:AbstractFloat, NAV<:AbstractNavigator} <: AbstractNavigatorState
    navigator::NAV                               # Navigator to be used
    topvol::Volume{T}                            # Typically the unplaced world
    currvol::Volume{T}                           # the current volume
    isinworld::Bool                              # inside world volume
    volstack::Vector{Int64}                      # keep the indexes of all daughters up to the current one 
    tolocal::Vector{Transformation3D{T}}         # Stack of transformations
end

#---Constructors------------------------------------------------------------------------------------
function NavigatorState{T}(top::Volume{T}, nav::AbstractNavigator=TrivialNavigator{T}(top)) where T
    x = NavigatorState{T,typeof(nav)}(nav, top, top, false, Vector{Int64}(), Vector{Transformation3D{T}}()) 
    sizehint!(x.volstack,16)
    sizehint!(x.tolocal,16)
    return x
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

#---Basic loops implemented using acceleration structures-------------------------------------------
@inline containedDaughters(::TrivialNavigator{T}, vol::Volume{T}, ::Point3{T}) where T = (i for i in 1:length(vol.daughters))
@inline function containedDaughters(nav::BVHNavigator{T}, vol::Volume{T}, point::Point3{T}) where T
    id = objectid(vol)
    if haskey(nav.bvhdict, id)
        bvh = nav.bvhdict[id]
        return (i for i in Iterators.flatten(inds for inds in pvolindices(x -> inside(point,x), bvh)))
    else
        return (i for i in 1:length(vol.daughters))
    end
end


#---Locate global point and initialize state-------------------------------------------------------- 
function locateGlobalPoint!(state::NavigatorState{T,NAV}, point::Point3{T}) where {T, NAV}
    nav = state.navigator
    reset!(state)
    vol = state.topvol
    if contains(vol, point)
        state.currvol = vol
        state.isinworld = true
    end
    # check the daughters in a recursive manner
    isinside = true
    while isinside
        isinside = false
        #for pvol in vol.daughters
        inds = containedDaughters(nav, vol, point)
        for i in inds
            pvol = vol.daughters[i]
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
