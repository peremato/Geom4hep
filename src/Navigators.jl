
abstract type AbstractNavigator end
abstract type AbstractNavigatorState end

#---Define the various navigators-------------------------------------------------------------------
struct TrivialNavigator{T<:AbstractFloat} <: AbstractNavigator
    world::Volume{T}   
end
function TrivialNavigator(world::Volume{T}) where T
    TrivialNavigator{T}(world)
end 

#---BVHNavigator-------------------------------------------------------------------------------------
struct BVHNavigator{T<:AbstractFloat} <: AbstractNavigator
    world::Volume{T}
    param::BVHParam{T}
    bvhdict::Dict{UInt64,BVH{T}}
    # Used to cache the lists of indices to placed volumes (avoid useless allocations)
    pvolind::Vector{Int64}
end

function BVHNavigator(world::Volume{T}; SAHfrac::T=T(8)) where T
    nav = BVHNavigator{T}(world, BVHParam{T}(SAHfrac), Dict{UInt64,BVH{T}}(),Int64[])
    #--- Create accelerationstructures---------------------------------
    buildBVH(nav, world, Set{UInt64}())
    sizehint!(nav.pvolind,1024)
    nav
end

function buildBVH(nav::BVHNavigator{T}, vol::Volume{T}, set::Set{UInt64}) where T
    id = vol.id
    #---return immediatelly if BVH is already there----
    id in set && return
    push!(set, id)
    #---only for more than 4 daughters
    if length(vol.daughters) > 4
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
    isinworld::Bool                              # inside world volume
    pvolstack::Stack{PlacedVolume{T}}            # stack of PlacedVolume
end

#---Constructors------------------------------------------------------------------------------------
function NavigatorState{T}(top::Volume{T}, nav::AbstractNavigator=TrivialNavigator{T}(top)) where T
    NavigatorState{T,typeof(nav)}(nav, top, false, Stack{PlacedVolume{T}}()) 
end

@inline function reset!(state::NavigatorState{T}) where T<:AbstractFloat
    empty!(state.pvolstack)
    state.isinworld = false 
end

@inline function currentVolume(state::NavigatorState{T}) where T<:AbstractFloat
    isempty(state.pvolstack) ? state.topvol : first(state.pvolstack).volume
end
@inline function currentTransformation(state::NavigatorState{T}) where T<:AbstractFloat
    isempty(state.pvolstack) ? one(Transformation3D{T}) : first(state.pvolstack).transformation
end

@inline function popOut!(state::NavigatorState{T}, point::Point3{T}) where T<:AbstractFloat
    if !isempty(state.pvolstack)
        point = point * currentTransformation(state) # transform to the mother system of reference
        pop!(state.pvolstack)
        relocatePoint!(state, point; inmother=true)
    else
        state.isinworld = false
    end
end

@inline function pushIn!(state::NavigatorState{T}, pvol::PlacedVolume{T}, point::Point3{T}) where T<:AbstractFloat
    push!(state.pvolstack, pvol)
    relocatePoint!(state, pvol.transformation * point)
end

#---Basic loops implemented using acceleration structures-------------------------------------------
containedDaughters(::TrivialNavigator{T}, vol::Volume{T}, ::Point3{T}) where T = eachindex(vol.daughters)
function containedDaughters(nav::BVHNavigator{T}, vol::Volume{T}, point::Point3{T}) where T
    id = vol.id
    if haskey(nav.bvhdict, id)
        bvh = nav.bvhdict[id]
        pvolindices!(nav.pvolind, x -> inside(x, point), bvh)
    else
        append!(empty!(nav.pvolind), 1:length(vol.daughters))
    end
    return nav.pvolind
end

intersectedDaughters(::TrivialNavigator{T}, vol::Volume{T}, ::Point3{T}, ::Vector3{T}) where T = eachindex(vol.daughters)
function intersectedDaughters(nav::BVHNavigator{T}, vol::Volume{T}, point::Point3{T}, dir::Vector3{T}) where T
    id = vol.id
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
        state.isinworld = true
    end
    # check the daughters in a recursive manner
    isinside = true
    while isinside
        isinside = false
        for i in containedDaughters(nav, vol, point)
            pvol = vol.daughters[i]
            if contains(pvol, point)
                push!(state.pvolstack, pvol)
                isinside = true
                vol   = pvol.volume
                point = pvol.transformation * point
                break
            end
        end
    end
    return state.isinworld
end


#---Locate global point and initialize state-------------------------------------------------------- 
function relocatePoint!(state::NavigatorState{T,NAV}, point::Point3{T}; inmother::Bool=false) where {T, NAV}
    nav  = state.navigator        # use the navigator acceleration
    vol  = currentVolume(state)   # start from data has been just updated
    # if last does not contain the point go up the hierarchy
    if inmother
        while !contains(vol, point)
            point = point * currentTransformation(state)
            if !isempty(state.pvolstack)
                pop!(state.pvolstack)
            end
            vol = currentVolume(state)
        end
    end
    # check the daughters in a non-recursive manner
    isinside = true
    while isinside
        isinside = false
        for i in containedDaughters(nav, vol, point)
            pvol = vol.daughters[i]
            if contains(pvol, point)
                push!(state.pvolstack, pvol)
                isinside = true
                vol   = pvol.volume
                point = pvol.transformation * point
                break
            end
        end
    end
    return nothing
end

@inline function isInVolume(state::NavigatorState{T}) where T<:AbstractFloat
    state.isinworld
end

function getClosestDaughter(nav::NAV, volume::Volume{T}, point::Point3{T}, dir::Vector3{T}, step_limit::T ) where {T,NAV}
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
    lpoint, ldir = transform((pvol.transformation for pvol in state.pvolstack.store), gpoint, gdir) 
    volume = currentVolume(state)
    step, idx = getClosestDaughter(state.navigator, volume, lpoint, ldir, step_limit)
    #---If didn't hit any daughter return distance to out
    if idx == 0
        step = distanceToOut(volume.shape, lpoint, ldir)
        if step > -kTolerance(T)
            popOut!(state, lpoint + ldir * (step + kPushTolerance(T)))  # go to mother or siblings
        elseif step < 0
            shape = volume.shape
            volstack = [pvol.idx for pvol in state.pvolstack]
            error("Negative distanceToOut. \nstep = $step \nshape = $shape \nvolstack = $volstack \npoint = $lpoint \ndir =   $ldir" )
        end
    else
    #---We hit a daughter, push it into the stack
        pvol = volume.daughters[idx]
        pushIn!(state, pvol, lpoint + ldir * (step + kPushTolerance(T))) # go to daughter or grand-daughters
    end
    return step
end
