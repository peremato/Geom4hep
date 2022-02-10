using Geom4hep
using CUDA
using BenchmarkTools, Test, Printf


@enum ShapeEnum::Int32 kBox kTube kTrd kCone
function getShapeType(shape::AbstractShape{T})::ShapeEnum where T<:AbstractFloat
    if typeof(shape) == Box{T}
        return kBox
    elseif typeof(shape) == Tube{T}
        return kTube
    elseif typeof(shape) == Trd{T}
        return kTrd
    elseif typeof(shape) == Cone{T}
        return kCone
    end
end 

struct CuVolume{T<:AbstractFloat}
    shapeTyp::ShapeEnum
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
    transformation::Transformation3D{T}
    volume::UInt32
end

struct CuModel{T<:AbstractFloat}
    indexes::Dict{UInt64, UInt32}
    volumes::Vector{CuVolume{T}}
    materials::Vector{CuMaterial{T}}
    daughters::Vector{CuPlacedVolume{T}}
    boxes::Vector{Box{T}}
    tubes::Vector{Tube{T}}
    trds::Vector{Trd{T}}
    cones::Vector{Cone{T}}
    function CuModel{T}() where T<:AbstractFloat
        new(Dict{UInt64, UInt32}(), 
            Vector{CuVolume{T}}(),
            Vector{CuMaterial{T}}(),
            Vector{CuPlacedVolume{T}}(),
            Vector{Box{T}}(),
            Vector{Tube{T}}(),
            Vector{Trd{T}}(),
            Vector{Cone{T}}())
    end
end

function pushObject(indexes::Dict{UInt64, UInt32}, vector::Vector{CUOBJ}, obj::OBJ)::UInt32 where {CUOBJ,OBJ}
    id = objectid(obj)
    if !haskey(indexes, id)
        push!(vector, CUOBJ == OBJ ? obj : CUOBJ(obj))
        indexes[id] = lastindex(vector)
    end
    indexes[id]
end

function pushVolume(model::CuModel{T}, vol::Volume{T}) where T <: AbstractFloat
    (; indexes, volumes, materials, daughters ) = model
    @show vol.label
    id = objectid(vol)
    @show haskey(indexes, id)
    if !haskey(indexes, id)
        volIdx = lastindex(model.volumes) + 1
        resize!(model.volumes, volIdx) 
        shapeTyp = getShapeType(vol.shape)
        shapeIdx = pushObject(indexes, getfield(model, 5 + Int(shapeTyp)), vol.shape)
        materialIdx = pushObject(indexes, materials, vol.material)
        daughterOff = lastindex(daughters)
        daughterLen = length(vol.daughters)
        @show daughterOff, daughterLen
        resize!(model.daughters, daughterOff + daughterLen)
        for d in 1:daughterLen
            model.daughters[d] = CuPlacedVolume{T}(vol.daughters[d].transformation, pushVolume(model, vol.daughters[d].volume ))
        end
        model.volumes[volIdx] = CuVolume{T}(shapeTyp, shapeIdx, materialIdx, daughterOff, daughterLen)    
        indexes[id] = volIdx
    end
    indexes[id]
end


if CUDA.functional()
    to_gpu(x::AbstractArray) = CuArray(x)
else
    to_gpu(x::AbstractArray) = x
end

end
function CuGeom(vol::Volume{T}) where T<:AbstractFloat
    model = CuModel{Float64}()
    pushVolume(model, vol)
    return (volumes = to_gpu(model.volumes), 
            materials = to_gpu(model.materials),
            daughters = to_gpu(model.daughters),
            boxes = to_gpu(model.boxes),
            tubes = to_gpu(model.tubes),
            trds = to_gpu(model.trds),
            cones = to_gpu(model.cones))         
end



function generateXRay(world::Volume{T}, npoints::Number, view::Symbol=:x) where T<:AbstractFloat
    lower, upper = extent(world.shape)
    dim = upper - lower
    if view == :x
        dim_a, dim_b = dim[2], dim[3]
    elseif view == :y
        dim_a, dim_b = dim[3], dim[1]
    elseif view == :z
        dim_a, dim_b = dim[1], dim[2]
    end
    pixel = round(sqrt(dim_a*dim_b/npoints), sigdigits=3)
    nx, ny = round(Int, dim_a/pixel), round(Int, dim_b/pixel)
    state = NavigatorState{Float64}(world)
    result = CUDA.zeros(nx,ny)
    threads = (16,16)
    blocks = cld.((nx,ny),threads)
    CUDA.@sync @cuda threads=threads blocks=blocks k_generateXRay(result, lower, pixel, state, view)
    return Array(result)
end


function k_generateXRay(result, lower::Point3{T}, pixel::T, state::NavigatorState{T}, view::Symbol) where T<:AbstractFloat
    nx, ny = size(result)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if i <= nx && j <= ny
        if view == :x
            point = Point3{T}(lower[1]+kTolerance(), lower[2]+(i-0.5)*pixel, lower[3]+(j-0.5)*pixel)
            dir = Vector3{T}(1,0,0)
        elseif view == :y
            point = Point3{T}(lower[1]+(j-0.5)*pixel, lower[2]+kTolerance(), lower[3]+(i-0.5)*pixel)
            dir = Vector3{T}(0,1,0)
        elseif shift == :z
            point = Point3{T}(lower[1]+(i-0.5)*pixel, lower[2]+(j-0.5)*pixel, lower[3]+kTolerance())
            dir = Vector3{T}(0,0,1)
        end
        locateGlobalPoint!(state, point)
        mass::T =  0.0
        step::T = -1.0
        while step != 0.0
            density = currentVolume(state).material.density
            step = computeStep!(state, point, dir, 1000.)
            point = point + dir * step
            mass += step * density
        end
        @inbounds result[i,j] = mass
    end
    return
end

#-----build and generate image-----------------------------
world = processGDML("examples/boxes.gdml")
