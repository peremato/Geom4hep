using LightXML
name = LightXML.name

include("Units.jl")

materials = Dict{String,AbstractMaterial}()
solids = Dict{String,AbstractShape{Float64}}()
volumes =  Dict{String,Volume{Float64}}()

# process <materials/>
function fillMaterials!(materials, element::XMLElement)
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            if elemname == "isotope"
                attrs = attributes_dict(e)
                matname = attrs["name"]
                mass = 0.
                for cc in child_nodes(e)
                    if name(cc) == "atom"
                        aa = attributes_dict(XMLElement(cc))
                        mass = parse(Float64, aa["value"])
                        break
                    end
                end
                materials[matname] = Isotope(matname, 
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
                        push!(composition,(fraction=parse(Float64, aa["n"]), isotope=materials[aa["ref"]]))
                    end
                end
                materials[matname] = Element(matname, composition)
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
                            args[:temperature] = parse(Float64, aa["value"])
                        elseif name(cc) == "D"
                            args[:density] = parse(Float64, aa["value"])
                        elseif name(cc) == "atom"
                            args[:mass] = parse(Float64, aa["value"])
                        elseif name(cc) == "fraction"
                            push!(composition,(fraction=parse(Float64, aa["n"]), element=materials[aa["ref"]]))
                            args[:composition] = composition
                        end
                    end
                end
                materials[matname] = Material(matname; args...)
            end
        end
    end
end

# process <solids/>
function fillSolids!(solids, element::XMLElement)
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "box"
                unit = eval(Meta.parse(attrs["lunit"]))
                solids[attrs["name"]] = Box{Float64}(parse(Float64, attrs["x"]) * unit, 
                                                     parse(Float64, attrs["y"]) * unit, 
                                                     parse(Float64, attrs["z"]) * unit)
            elseif elemname == "trd"
                unit = eval(Meta.parse(attrs["lunit"]))
                solids[attrs["name"]] = Trd{Float64}(parse(Float64, attrs["x1"]) * unit, 
                                                     parse(Float64, attrs["x2"]) * unit, 
                                                     parse(Float64, attrs["y1"]) * unit,
                                                     parse(Float64, attrs["y2"]) * unit,
                                                     parse(Float64, attrs["z"]) * unit)
            else
                solids[attrs["name"]] = NoShape{Float64}(elemname)
            end
        end
    end
end

# process <structure/>
function fillVolumes!(volumes, element::XMLElement)
    for c in child_nodes(element) #--- loop over volumes
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "volume"
                volname = attrs["name"]
                shape = nothing
                material = nothing
                volume = nothing
                for cc in child_nodes(e)  #--- loop over ref to solid, material and daughters
                    if is_elementnode(cc)
                        aa = attributes_dict(XMLElement(cc))
                        if  name(cc) == "solidref"
                            shape = solids[aa["ref"]]
                            if !isnothing(material)
                                volume = Volume{Float64}(volname, shape, material)
                            end
                        elseif name(cc) == "materialref"
                            material = materials[aa["ref"]]
                            if !isnothing(shape)
                                volume = Volume{Float64}(volname, shape, material)
                            end
                        elseif name(cc) == "physvol"
                            daughter = nothing
                            position = (0.,0.,0.)
                            rotation = (0.,0.,0.)
                            for ccc in child_nodes(XMLElement(cc))
                                if is_elementnode(ccc)
                                    aa = attributes_dict(XMLElement(ccc))
                                    if name(ccc) == "volumeref"
                                        daughter = volumes[aa["ref"]]
                                    elseif name(ccc) == "position"
                                        unit = eval(Meta.parse(aa["unit"]))
                                        position = (parse(Float64, aa["x"]) * unit, 
                                                    parse(Float64, aa["y"]) * unit,
                                                    parse(Float64, aa["z"]) * unit)
                                    elseif name(ccc) == "rotation"
                                        unit = eval(Meta.parse(aa["unit"]))
                                        rotation = (parse(Float64, aa["x"]) * unit, 
                                                    parse(Float64, aa["y"]) * unit,
                                                    parse(Float64, aa["z"]) * unit)
                                    end
                                end
                            end
                            placeDaughter!(volume, Transformation3D{Float64}(position..., rotation...), daughter)
                        end
                    end
                end
                volumes[volname] = volume
            end
        end
    end
end


# process the full file
function processGDML(fname::String)
    world = nothing
    xdoc = parse_file(fname)
    xroot = root(xdoc) 
    if name(xroot) != "gdml"
        throw(ErrorException("File @fname is not a GDML file"))
    end
    empty!(materials)
    empty!(volumes)
    empty!(solids)
    for c in child_nodes(xroot)  # c is an instance of XMLNode
        if is_elementnode(c)
            e = XMLElement(c)    # this makes an XMLElement instance
            if name(e) == "materials"
                fillMaterials!(materials, e)
            elseif name(e) == "solids"
                fillSolids!(solids, e)
            elseif name(e) == "structure"
                fillVolumes!(volumes, e)
            elseif name(e) == "setup"
                for cc in child_nodes(e)
                    if is_elementnode(cc)
                        ee = XMLElement(cc)
                        attrs = attributes_dict(ee)
                        if name(ee) == "world"
                           world = volumes[attrs["ref"]] 
                        end
                    end
                end
            end
        end
    end
    free(xdoc)
    return world
end

