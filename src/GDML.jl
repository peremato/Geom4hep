using LightXML, Printf
name = LightXML.name

include("Units.jl")

struct GDMLDicts{T<:AbstractFloat, S}
    materials::Dict{String,AbstractMaterial{T}}
    solids::Dict{String,AbstractShape{T}}
    volumes::Dict{String,VolumeP{T, S}}
    positions::Dict{String, Vector3{T}}
    rotations::Dict{String, RotXYZ{T}}
    function GDMLDicts{T, S}() where {T<:AbstractFloat, S}
        new(Dict{String,AbstractMaterial{T}}(),
            Dict{String,AbstractShape{T}}(),
            Dict{String,VolumeP{T, S}}(),
            Dict{String, Vector3{T}}(),
            Dict{String, RotXYZ{T}}()
        )
    end
end

# process <define/>
function fillDefine!(dicts::GDMLDicts{T, S}, element::XMLElement) where {T<:AbstractFloat, S}
    (; positions, rotations) = dicts
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "position"
                unit = eval(Meta.parse(attrs["unit"]))
                positions[attrs["name"]] = Vector3{T}(parse(T, attrs["x"]) * unit, 
                                                      parse(T, attrs["y"]) * unit, 
                                                      parse(T, attrs["z"]) * unit)
            elseif elemname == "rotation"
                unit = eval(Meta.parse(attrs["unit"]))
                rotations[attrs["name"]] = RotZYX{T}(parse(T, attrs["z"]) * unit, 
                                                     parse(T, attrs["y"]) * unit, 
                                                     parse(T, attrs["x"]) * unit)
            else
                @printf "define %s not yet suported\n" elemname
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
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "box"
                lunit = eval(Meta.parse(attrs["lunit"]))
                solids[attrs["name"]] = Box{T}(parse(T, attrs["x"]) * lunit / 2, 
                                               parse(T, attrs["y"]) * lunit / 2, 
                                               parse(T, attrs["z"]) * lunit / 2)
            elseif elemname == "trd"
                lunit = eval(Meta.parse(attrs["lunit"]))
                solids[attrs["name"]] = Trd{T}(parse(T, attrs["x1"]) * lunit / 2, 
                                               parse(T, attrs["x2"]) * lunit / 2, 
                                               parse(T, attrs["y1"]) * lunit / 2,
                                               parse(T, attrs["y2"]) * lunit / 2,
                                               parse(T, attrs["z"])  * lunit / 2)
            elseif elemname == "trap"
                lunit = eval(Meta.parse(attrs["lunit"]))
                aunit = eval(Meta.parse(attrs["aunit"]))
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
                lunit = eval(Meta.parse(attrs["lunit"]))
                aunit = eval(Meta.parse(attrs["aunit"]))
                solids[attrs["name"]] = Tube{T}(parse(T, attrs["rmin"]) * lunit, 
                                                parse(T, attrs["rmax"]) * lunit, 
                                                parse(T, attrs["z"]) * lunit / 2,
                                                parse(T, attrs["startphi"]) * aunit,
                                                parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "cone"
                lunit = eval(Meta.parse(attrs["lunit"]))
                aunit = eval(Meta.parse(attrs["aunit"]))
                solids[attrs["name"]] = Cone{T}(parse(T, attrs["rmin1"]) * lunit, 
                                                parse(T, attrs["rmax1"]) * lunit,
                                                parse(T, attrs["rmin2"]) * lunit,
                                                parse(T, attrs["rmax2"]) * lunit, 
                                                parse(T, attrs["z"]) * lunit / 2,
                                                parse(T, attrs["startphi"]) * aunit,
                                                parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "polycone"
                lunit = eval(Meta.parse(attrs["lunit"]))
                aunit = eval(Meta.parse(attrs["aunit"]))
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
                lunit = eval(Meta.parse(attrs["lunit"]))
                aunit = eval(Meta.parse(attrs["aunit"]))
                solids[attrs["name"]] = CutTube{T}(parse(T, attrs["rmin"]) * lunit,
                                                   parse(T, attrs["rmax"]) * lunit,
                                                   parse(T, attrs["z"]) * lunit,
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
                            first = solids[aa["ref"]]
                        elseif name(cc) == "second"
                            second = solids[aa["ref"]]
                        end
                    end
                end
                if elemname == "union"
                    solids[attrs["name"]] = BooleanUnion(first, second, transformation)
                elseif elemname == "subtraction"
                    solids[attrs["name"]] = BooleanSubtraction(first, second, transformation)
                elseif elemname == "intersection"
                    solids[attrs["name"]] = BooleanIntersection(first, second, transformation)
                end
            else
                @printf "Shape %s not yet suported. Using NoShape\n" elemname
                solids[attrs["name"]] = NoShape{T}()
            end
        end
    end
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
                volume = elemname == "volume" ? VolumeP(volname, shape, material) : Assembly{T}(volname)
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
                unit = eval(Meta.parse(aa["unit"]))
                position = Vector3{T}(parse(T, aa["x"]) * unit, 
                                      parse(T, aa["y"]) * unit,
                                      parse(T, aa["z"]) * unit)
            elseif name(c) == "rotation"
                unit = eval(Meta.parse(aa["unit"]))
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
    dicts = GDMLDicts{T, PlacedVolume{T}}()
    world = nothing
    xdoc = parse_file(fname)
    xroot = root(xdoc) 
    if name(xroot) != "gdml"
        throw(ErrorException("File @fname is not a GDML file"))
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

