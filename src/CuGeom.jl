using CUDA, Adapt
export fillCuGeometry, CuNavigatorState, CuGeoModel

struct CuVolume{T<:AbstractFloat}
    shapeIdx::UInt32
    materialIdx::UInt32
    daughterOff::UInt32
    daughterLen::UInt32
end

struct CuMaterial{T<:AbstractFloat}
    density::T
    temperature::T
    Amass::T
end

CuMaterial{T}(m::Material) where T<: AbstractFloat = CuMaterial{T}(m.density, m.temperature, m.Amass) 

struct CuPlacedVolume{T<:AbstractFloat}
    idx::UInt32
    transformation::Transformation3D{T}
    volume::UInt32
end

struct CuGeoModel{V, M, PV, S}
    volumes::V
    materials::M
    placedvolumes::PV
    shapes::S
end

Adapt.@adapt_structure CuGeoModel

function CuGeoModel{T}() where T<:AbstractFloat
    CuGeoModel(
        Vector{CuVolume{T}}(),
        Vector{CuMaterial{T}}(),
        Vector{CuPlacedVolume{T}}(),
        Vector{Union{Box{T},Tube{T},Trd{T}}}()
    )
end

function Base.empty!(model::CuGeoModel) where T<:AbstractFloat
    for n in fieldnames(typeof(model))
        empty!(getfield(model, n))
    end
end

getShape(model::CuGeoModel, vol::CuVolume{T}) where T<:AbstractFloat = model.shapes[vol.shapeIdx]
contains(model::CuGeoModel, vol::CuVolume{T}, point::Point3{T}) where T<:AbstractFloat = inside(getShape(model, vol), point) == kInside
contains(model::CuGeoModel, pvol::CuPlacedVolume, point::Point3{T}) where T<:AbstractFloat = inside(getShape(model, model.volumes[pvol.volume]), pvol.transformation * point) == kInside
distanceToIn(model::CuGeoModel, pvol::CuPlacedVolume{T}, p::Point3{T}, d::Vector3{T}) where T<:AbstractFloat = distanceToIn(getShape(model, model.volumes[pvol.volume]), pvol.transformation * p, pvol.transformation * d)

mutable struct CuNavigatorState{T<:AbstractFloat}
    topVolume::UInt32
    currentVol::UInt32
    currentDepth::Int64
    volstack::SVector{16,UInt32}
    tolocal::Transformation3D{T}
    function CuNavigatorState{T}(top::Integer) where T<:AbstractFloat
        state = new{T}(top)
        reset!(state)
        return state
    end
end

@inline function reset!(state::CuNavigatorState{T}) where T<:AbstractFloat
    state.currentVol = state.topVolume
    state.currentDepth = 0
    state.volstack = zero(SVector{16,UInt32})
    state.tolocal = one(Transformation3D{T})
end

@inline function pushIn!(state::CuNavigatorState{T}, pvol::CuPlacedVolume{T}) where T<:AbstractFloat
    state.currentVol = pvol.volume
    state.currentDepth += 1
    state.volstack = setindex(state.volstack, pvol.idx, state.currentDepth)
    state.tolocal *= pvol.transformation
    nothing
end

@inline function popOut!(state::CuNavigatorState{T}, model::CuGeoModel) where T<:AbstractFloat
    if state.currentDepth > 0
        state.currentDepth -= 1
        #---Start from the begining to get the current volume (mother) and transformation
        state.currentVol = state.topVolume
        state.tolocal = one(Transformation3D{T})
        for depth in 1:state.currentDepth
            vol  = model.volumes[state.currentVol]
            pvol = model.placedvolumes[vol.daughterOff +  state.volstack[depth]]
            state.currentVol = pvol.volume
            state.tolocal *= pvol.transformation
        end
    end
    nothing
end

@inline function collectDaughters!(model::CuGeoModel, state::CuNavigatorState{T}, vol::CuVolume, point::Point3{T}) where T<:AbstractFloat
    isinside = true
    while isinside
        isinside = false
        for d in 1:vol.daughterLen
            pvol = model.placedvolumes[vol.daughterOff + d]
            if contains(model, pvol, point)
                pushIn!(state, pvol)
                # emulation of recursion
                # collectDaughters!(model, state, model.volumes[pvol.volume], pvol.transformation * point)
                isinside = true
                vol   = model.volumes[pvol.volume]
                point = pvol.transformation * point
                break
            end
        end
    end
end

@inline function locateGlobalPoint!(model::CuGeoModel, state::CuNavigatorState{T}, gpoint::Point3{T}) where T<:AbstractFloat
    state.currentVol = state.topVolume
    state.currentDepth = 0
    state.tolocal = one(Transformation3D{T})
    vol = model.volumes[state.currentVol]
    if contains(model, vol, gpoint)
        collectDaughters!(model, state, vol, gpoint)
        return true
    end
    return false
end

function getClosestDaughter(model::CuGeoModel, vol::CuVolume{T}, point::Point3{T}, dir::Vector3{T}, step_limit::T ) where T<:AbstractFloat
    step = step_limit
    candidate = 0
    #---Linear loop over all daughters
    for d in 1:vol.daughterLen
        pvol = model.placedvolumes[vol.daughterOff + d]
        #---Assuming that it is not yet inside the daughter (otherwise it returns -1.)
        dist = distanceToIn(model, pvol, point, dir)
        if dist > 0. &&dist != Inf && dist < step 
            step = dist
            candidate = d
        end
    end
    return step, candidate
end


@inline function computeStep!(model::CuGeoModel, state::CuNavigatorState{T}, gpoint::Point3{T}, gdir::Vector3{T}, step_limit::T) where T<:AbstractFloat
    lpoint, ldir = state.tolocal * gpoint, state.tolocal * gdir
    volume = model.volumes[state.currentVol]

    step, idx = getClosestDaughter(model, volume, lpoint, ldir, step_limit)

    #---If didn't hit any daughter return distance to out
    if idx == 0
        step = distanceToOut(getShape(model, volume), lpoint, ldir)
        if step > 0.
            popOut!(state, model)
        elseif step == 0.
            if state.currentDepth > 0
                popOut!(state, model)
            else
                step = -1.
            end
        end
    end
    #---We hit a daughter, push it into the stack
    if idx != 0
        step += kTolerance(T) # to ensure that do not stay in the surface of the daughter
        pvol = model.placedvolumes[volume.daughterOff + idx]
        pushIn!(state, pvol)
    end
    return step
end

#-----------------------------------------------------------------------------------------------------------------------
#---Transform the geometry model to Cuda model
#-----------------------------------------------------------------------------------------------------------------------
function pushObject(indexes::Dict{UInt64, UInt32}, vector::Vector{CUOBJ}, obj::OBJ)::UInt32 where {CUOBJ,OBJ}
    id = objectid(obj)
    if !haskey(indexes, id)
        push!(vector, CUOBJ == OBJ || typeof(CUOBJ) == Union ? obj : CUOBJ(obj))
        indexes[id] = lastindex(vector)
    end
    indexes[id]
end

function pushVolume!(indexes::Dict{UInt64, UInt32}, model::CuGeoModel, vol::Volume{T}) where T <: AbstractFloat
    (; volumes, materials, placedvolumes, shapes ) = model
    id = objectid(vol)
    if !haskey(indexes, id)
        volIdx = lastindex(model.volumes) + 1
        resize!(model.volumes, volIdx) 
        shapeIdx = pushObject(indexes, shapes, vol.shape)
        materialIdx = pushObject(indexes, materials, vol.material)
        daughterOff = lastindex(placedvolumes)
        daughterLen = length(vol.daughters)
        resize!(model.placedvolumes, daughterOff + daughterLen)
        for d in 1:daughterLen
            model.placedvolumes[d + daughterOff] = CuPlacedVolume{T}(d, vol.daughters[d].transformation, pushVolume!(indexes, model, vol.daughters[d].volume ))
        end
        model.volumes[volIdx] = CuVolume{T}(shapeIdx, materialIdx, daughterOff, daughterLen)    
        indexes[id] = volIdx
    end
    indexes[id]
end

if CUDA.functional()
    to_gpu(x::AbstractArray) = CuArray(x)
else
    to_gpu(x::AbstractArray) = x
end

function fillCuGeometry(vol::Volume{T}) where T<:AbstractFloat
    indexes = Dict{UInt64, UInt32}()
    model = CuGeoModel{T}()
    pushVolume!(indexes, model, vol)
    return model
end
