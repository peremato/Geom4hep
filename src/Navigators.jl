
abstract type AbstractNavigator end
abstract type AbstractNavigatorState end

#---Define the various navigators-------------------------------------------------------------------
struct TrivialNavigator{T<:AbstractFloat, S} <: AbstractNavigator
    world::VolumeP{T, S}
end
function TrivialNavigator(world::VolumeP{T, PV}) where {T, PV}
    TrivialNavigator{T}(world)
end 

#---BVHNavigator-------------------------------------------------------------------------------------
struct BVHNavigator{T<:AbstractFloat, S} <: AbstractNavigator
    world::VolumeP{T, S}
    param::BVHParam{T}
    bvhdict::Dict{UInt64,BVH{T}}
    # Used to cache the lists of indices to placed volumes (avoid useless allocations)
    pvolind::Vector{Int64}
end

function BVHNavigator(world::VolumeP{T, S}; SAHfrac::T=T(8)) where {T, S}
    nav = BVHNavigator{T, S}(world, BVHParam{T}(SAHfrac), Dict{UInt64,BVH{T}}(),Int64[])
    #--- Create accelerationstructures---------------------------------
    buildBVH(nav, world, Set{UInt64}())
    sizehint!(nav.pvolind,1024)
    nav
end

function buildBVH(nav::BVHNavigator{T}, vol::VolumeP{T, S}, set::Set{UInt64}) where {T, S}
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
mutable struct NavigatorState{T<:AbstractFloat, S, NAV<:AbstractNavigator} <: AbstractNavigatorState
    navigator::NAV                               # Navigator to be used
    topvol::VolumeP{T,S}                            # Typically the unplaced world
    currvol::VolumeP{T,S}                           # the current volume
    isinworld::Bool                              # inside world volume
    volstack::Vector{Int64}                      # keep the indexes of all daughters up to the current one 
    tolocal::Vector{Transformation3D{T}}         # Stack of transformations
end

#---Constructors------------------------------------------------------------------------------------
function NavigatorState{T}(top::VolumeP{T, S}, nav::AbstractNavigator=TrivialNavigator{T}(top)) where {T, S}
    x = NavigatorState(nav, top, top, false, Vector{Int64}(), Vector{Transformation3D{T}}()) 
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
containedDaughters(::TrivialNavigator{T}, vol::VolumeP, ::Point3{T}) where T = eachindex(vol.daughters)
function containedDaughters(nav::BVHNavigator{T}, vol::VolumeP, point::Point3{T}) where T
    id = objectid(vol)
    if haskey(nav.bvhdict, id)
        bvh = nav.bvhdict[id]
        pvolindices!(nav.pvolind, x -> inside(x, point), bvh)
    else
        append!(empty!(nav.pvolind), 1:length(vol.daughters))
    end
    return nav.pvolind
end

intersectedDaughters(::TrivialNavigator{T}, vol::VolumeP, ::Point3{T}, ::Vector3{T}) where T = eachindex(vol.daughters)
function intersectedDaughters(nav::BVHNavigator{T}, vol::VolumeP, point::Point3{T}, dir::Vector3{T}) where T
    id = objectid(vol)
    if haskey(nav.bvhdict, id)
        bvh = nav.bvhdict[id]
        pvolindices!(nav.pvolind, x -> intersect(x, point, dir), bvh)
    else
        append!(empty!(nav.pvolind), 1:length(vol.daughters))
    end
    return nav.pvolind
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
        for i in containedDaughters(nav, vol, point)
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

function getClosestDaughter(nav::NAV, volume::VolumeP, point::Point3{T}, dir::Vector3{T}, step_limit::T ) where {T,NAV}
    step = step_limit
    candidate = 0
    #---Linear loop over the daughters-------------------------------------------------
    for i in intersectedDaughters(nav, volume, point, dir)
        pvol = volume.daughters[i]
        #---Assuming that it is not yet inside the daughter (otherwise it returns -1.)
        dist = distanceToIn(pvol, point, dir)
        dist < 0. && return (zero(T), pvol.idx ) 
        if dist > 0. && dist != Inf && dist < step
            step = dist
            candidate = pvol.idx
        end
    end
    return step, candidate
end

#---Update the NavigatorState and get the step (distance)
function computeStep!(state::NavigatorState{T,NAV}, gpoint::Point3{T}, gdir::Vector3{T}, step_limit::T) where {T,NAV}
    lpoint, ldir = transform(state.tolocal, gpoint, gdir) 
    volume = currentVolume(state)
    step, idx = getClosestDaughter(state.navigator, volume, lpoint, ldir, step_limit)
    #---If didn't hit any daughter return distance to out
    if idx == 0
        step = distanceToOut(volume.shape, lpoint, ldir)
        if step > -kTolerance(T)
            popOut!(state)
        elseif step < 0
            shape = volume.shape
            volstack = state.volstack
            error("Negative distanceToOut. \nstep = $step \nshape = $shape \nvolstack = $volstack \npoint = $lpoint \ndir =   $ldir" )
        end
    else
    #---We hit a daughter, push it into the stack
        step += kPushTolerance(T)  # to ensure that do not stay in the surface of the daughter
        pvol = volume.daughters[idx]
        pushIn!(state, pvol)
    end
    return step
end
