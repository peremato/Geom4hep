using LightXML, Printf
const name = LightXML.name

#--- Unit lookup used when reading GDML (avoids eval of attribute strings)
using .Units: mm, cm, m, km, um, nm, rad, radian, mrad, deg, degree
const _length_units = Dict{String,Float64}("mm" => mm, "cm" => cm, "m" => m, "km" => km,
                                           "um" => um, "micrometer" => um, "nm" => nm)
const _angle_units  = Dict{String,Float64}("rad" => rad, "radian" => radian, "mrad" => mrad,
                                           "deg" => deg, "degree" => degree)

function _parse_unit(table::Dict{String,Float64}, unit::AbstractString)
    haskey(table, unit) || throw(ArgumentError("Unknown or unsupported unit \"$unit\""))
    table[unit]
end

"""
    _unit(attrs, key, default)

Unit factor for the GDML attribute `key` (e.g. "lunit", "aunit", "unit"). Missing attributes take
the GDML defaults (mm for lengths, rad for angles).
"""
function _unit(attrs::AbstractDict, key::AbstractString, default::AbstractString)
    u = get(attrs, key, default)
    haskey(_length_units, u) && return _length_units[u]
    _parse_unit(_angle_units, u)
end
_lunit(attrs) = _unit(attrs, "lunit", "mm")
_aunit(attrs) = _unit(attrs, "aunit", "rad")

struct GDMLDicts{T<:AbstractFloat}
    materials::Dict{String,AbstractMaterial{T}}
    solids::Dict{String,AbstractShape{T}}
    volumes::Dict{String,Volume{T}}
    positions::Dict{String, Vector3{T}}
    rotations::Dict{String, RotXYZ{T}}
    function GDMLDicts{T}() where T<:AbstractFloat
        new(Dict{String,AbstractMaterial{T}}(),
            Dict{String,AbstractShape{T}}(),
            Dict{String,Volume{T}}(),
            Dict{String, Vector3{T}}(),
            Dict{String, RotXYZ{T}}()
        )
    end
end

# process <define/>
function fillDefine!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; positions, rotations) = dicts
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "position"
                unit = _unit(attrs, "unit", "mm")
                positions[attrs["name"]] = Vector3{T}(parse(T, attrs["x"]) * unit, 
                                                      parse(T, attrs["y"]) * unit, 
                                                      parse(T, attrs["z"]) * unit)
            elseif elemname == "rotation"
                unit = _unit(attrs, "unit", "rad")
                rotations[attrs["name"]] = RotZYX{T}(parse(T, attrs["z"]) * unit, 
                                                     parse(T, attrs["y"]) * unit, 
                                                     parse(T, attrs["x"]) * unit)
            else
                @warn "GDML: define $elemname not yet supported" maxlog=1
            end
        end
    end
end

# process <materials/>
function fillMaterials!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; materials) = dicts
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            if elemname == "isotope"
                attrs = attributes_dict(e)
                matname = attrs["name"]
                mass = zero(T)
                for cc in child_nodes(e)
                    if name(cc) == "atom"
                        aa = attributes_dict(XMLElement(cc))
                        mass = parse(T, aa["value"])
                        break
                    end
                end
                materials[matname] = Isotope{T}(matname, 
                                                parse(Int32, attrs["N"]), 
                                                parse(Int32, attrs["Z"]), 
                                                mass)
            elseif elemname == "element"
                attrs = attributes_dict(e)
                matname = attrs["name"]
                composition = []
                for cc in child_nodes(e)
                    if name(cc) == "fraction"
                        aa = attributes_dict(XMLElement(cc))
                        push!(composition,(fraction=parse(T, aa["n"]), isotope=materials[aa["ref"]]))
                    end
                end
                materials[matname] = Element{T}(matname, composition)
            elseif elemname == "material"
                attrs = attributes_dict(e)
                matname = attrs["name"]
                args = Dict{Symbol, Any}()
                composition = []
                if haskey(attrs,"state")
                    args[:state] = attrs["state"]
                end
                for cc in child_nodes(e)
                    if is_elementnode(cc)
                        aa = attributes_dict(XMLElement(cc))
                        if name(cc) == "T"
                            args[:temperature] = parse(T, aa["value"])
                        elseif name(cc) == "D"
                            args[:density] = parse(T, aa["value"])
                        elseif name(cc) == "atom"
                            args[:mass] = parse(T, aa["value"])
                        elseif name(cc) == "fraction"
                            push!(composition,(fraction=parse(T, aa["n"]), element=materials[aa["ref"]]))
                            args[:composition] = composition
                        end
                    end
                end
                materials[matname] = Material{T}(matname; args...)
            end
        end
    end
end

# process <solids/>
function fillSolids!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; solids) = dicts
    unsupported = Dict{String,Int}()
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "box"
                lunit = _lunit(attrs)
                solids[attrs["name"]] = Box{T}(parse(T, attrs["x"]) * lunit / 2, 
                                               parse(T, attrs["y"]) * lunit / 2, 
                                               parse(T, attrs["z"]) * lunit / 2)
            elseif elemname == "trd"
                lunit = _lunit(attrs)
                solids[attrs["name"]] = Trd{T}(parse(T, attrs["x1"]) * lunit / 2, 
                                               parse(T, attrs["x2"]) * lunit / 2, 
                                               parse(T, attrs["y1"]) * lunit / 2,
                                               parse(T, attrs["y2"]) * lunit / 2,
                                               parse(T, attrs["z"])  * lunit / 2)
            elseif elemname == "trap"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = Trap{T}(parse(T, attrs["z"]) * lunit / 2, 
                                                parse(T, attrs["theta"]) * aunit, 
                                                parse(T, attrs["phi"]) * aunit,
                                                parse(T, attrs["y1"]) * lunit / 2,
                                                parse(T, attrs["x1"])  * lunit / 2,
                                                parse(T, attrs["x2"])  * lunit / 2,
                                                parse(T, attrs["alpha1"])  * aunit,
                                                parse(T, attrs["y2"])  * lunit / 2,
                                                parse(T, attrs["x3"])  * lunit / 2,
                                                parse(T, attrs["x4"])  * lunit / 2,
                                                parse(T, attrs["alpha2"])  * aunit)
            elseif elemname == "tube"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = Tube{T}(parse(T, attrs["rmin"]) * lunit, 
                                                parse(T, attrs["rmax"]) * lunit, 
                                                parse(T, attrs["z"]) * lunit / 2,
                                                parse(T, attrs["startphi"]) * aunit,
                                                parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "cone"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = Cone{T}(parse(T, attrs["rmin1"]) * lunit, 
                                                parse(T, attrs["rmax1"]) * lunit,
                                                parse(T, attrs["rmin2"]) * lunit,
                                                parse(T, attrs["rmax2"]) * lunit, 
                                                parse(T, attrs["z"]) * lunit / 2,
                                                parse(T, attrs["startphi"]) * aunit,
                                                parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "polycone"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                rmax, rmin, z = Vector{T}(), Vector{T}(), Vector{T}()
                for cc in child_nodes(e)
                    if name(cc) == "zplane"
                        aa = attributes_dict(XMLElement(cc))
                        push!(rmax, parse(T, aa["rmax"]) * lunit)
                        push!(rmin, parse(T, aa["rmin"]) * lunit)
                        push!(z, parse(T, aa["z"]) * lunit)
                    end
                end
                N = length(rmax)
                solids[attrs["name"]] = Polycone{T}(rmin, rmax, z,
                                                    parse(T, attrs["startphi"]) * aunit,
                                                    parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "cutTube"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = CutTube{T}(parse(T, attrs["rmin"]) * lunit,
                                                   parse(T, attrs["rmax"]) * lunit,
                                                   parse(T, attrs["z"]) * lunit / 2,   # GDML gives the full length
                                                   parse(T, attrs["startphi"]) * aunit,
                                                   parse(T, attrs["deltaphi"]) * aunit,
                                                   Vector3{T}(parse(T, attrs["lowX"]), parse(T, attrs["lowY"]), parse(T, attrs["lowZ"])),
                                                   Vector3{T}(parse(T, attrs["highX"]), parse(T, attrs["highY"]), parse(T, attrs["highZ"])))
            elseif elemname in ("union", "subtraction", "intersection") 
                first = nothing
                second = nothing
                transformation = getTransformation(dicts, e)
                for cc in child_nodes(e)
                    if is_elementnode(cc)
                        aa = attributes_dict(XMLElement(cc))
                        if name(cc) == "first"
                            first = _get_solid(solids, aa["ref"], attrs["name"])
                        elseif name(cc) == "second"
                            second = _get_solid(solids, aa["ref"], attrs["name"])
                        end
                    end
                end
                if isnothing(first) || isnothing(second)
                    # Broken GDML: keep whatever operand is valid (or nothing at all for intersections)
                    solids[attrs["name"]] = elemname == "intersection" || isnothing(first) ? NoShape{T}() : first
                elseif elemname == "union"
                    solids[attrs["name"]] = BooleanUnion(first, second, transformation)
                elseif elemname == "subtraction"
                    solids[attrs["name"]] = BooleanSubtraction(first, second, transformation)
                elseif elemname == "intersection"
                    solids[attrs["name"]] = BooleanIntersection(first, second, transformation)
                end
            else
                unsupported[elemname] = get(unsupported, elemname, 0) + 1
                solids[attrs["name"]] = NoShape{T}()
            end
        end
    end
    isempty(unsupported) || @warn "GDML: unsupported shapes replaced by NoShape" unsupported
end

function _get_solid(solids::Dict, ref::AbstractString, owner::AbstractString)
    haskey(solids, ref) && return solids[ref]
    @warn "GDML: solid \"$ref\" referenced by boolean \"$owner\" is not defined; ignoring this operand"
    return nothing
end

# process <structure/>
function fillVolumes!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; materials, solids, volumes, positions, rotations) = dicts
    for c in child_nodes(element) #--- loop over volumes
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname in ("volume", "assembly")
                volname = attrs["name"]
                shape = nothing
                material = nothing
                pvols = PlacedVolume{T}[]
                for cc in child_nodes(e)  #--- loop over ref to solid, material and daughters
                    if is_elementnode(cc)
                        aa = attributes_dict(XMLElement(cc))
                        if  name(cc) == "solidref"
                            shape = solids[aa["ref"]]
                        elseif name(cc) == "materialref"
                            material = materials[aa["ref"]]
                        elseif name(cc) == "physvol"
                            daughter = nothing
                            transformation = getTransformation(dicts, XMLElement(cc))
                            for ccc in child_nodes(XMLElement(cc))
                                if is_elementnode(ccc)
                                    aa = attributes_dict(XMLElement(ccc))
                                    if name(ccc) == "volumeref"
                                        daughter = volumes[aa["ref"]]
                                    end
                                end
                            end
                            push!(pvols, PlacedVolume{T}(-1, transformation, daughter))
                        end
                    end
                end
                volume = elemname == "volume" ? Volume{T}(volname, shape, material) : Assembly{T}(volname)
                for pvol in pvols
                    placeDaughter!(volume, pvol.transformation, pvol.volume)
                end
                volumes[volname] = volume
            end
        end
    end
end

# fill Transformation3D
function getTransformation(dicts::GDMLDicts{T}, element::XMLElement) where T
    (; positions, rotations) = dicts
    position = Vector3{T}(0,0,0)
    rotation = RotZYX{T}(0,0,0)
    for c in child_nodes(element)
        if is_elementnode(c)
            aa = attributes_dict(XMLElement(c))
            if name(c) == "position"
                unit = _unit(aa, "unit", "mm")
                position = Vector3{T}(parse(T, aa["x"]) * unit, 
                                      parse(T, aa["y"]) * unit,
                                      parse(T, aa["z"]) * unit)
            elseif name(c) == "rotation"
                unit = _unit(aa, "unit", "rad")
                rotation = RotZYX{T}(parse(T, aa["z"]) * unit, 
                                     parse(T, aa["y"]) * unit, 
                                     parse(T, aa["x"]) * unit)
            elseif name(c) == "positionref"
                position = positions[aa["ref"]]
            elseif name(c) == "rotationref"
                rotation = rotations[aa["ref"]]
            end
        end
    end
    Transformation3D{T}(RotMatrix3{T}(rotation), position)
end


# process the full file
function processGDML(fname::String, ::Type{T}=Float64) where T
    dicts = GDMLDicts{T}()
    world = nothing
    xdoc = parse_file(fname)
    xroot = root(xdoc) 
    if name(xroot) != "gdml"
        throw(ArgumentError("File $fname is not a GDML file"))
    end
    for c in child_nodes(xroot)  # c is an instance of XMLNode
        if is_elementnode(c)
            e = XMLElement(c)    # this makes an XMLElement instance
            if name(e) == "materials"
                fillMaterials!(dicts, e)
            elseif name(e) == "solids"
                fillSolids!(dicts, e)
            elseif name(e) == "structure"
                fillVolumes!(dicts, e)
            elseif name(e) == "define"
                fillDefine!(dicts, e)
            elseif name(e) == "setup"
                for cc in child_nodes(e)
                    if is_elementnode(cc)
                        ee = XMLElement(cc)
                        attrs = attributes_dict(ee)
                        if name(ee) == "world"
                           world = dicts.volumes[attrs["ref"]] 
                        end
                    end
                end
            end
        end
    end
    free(xdoc)
    return world
end

