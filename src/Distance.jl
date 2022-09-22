function distanceToIn(shape, point, dir)
    if shape isa Trap
        distanceToIn_trap(shape, point, dir)
    elseif shape isa Trd
        distanceToIn_trd(shape, point, dir)
    elseif shape isa Cone
        distanceToIn_cone(shape, point, dir)
    elseif shape isa Box
        distanceToIn_box(shape, point, dir)
    elseif shape isa Tube
        distanceToIn_tube(shape, point, dir)
    elseif shape isa Polycone
        distanceToIn_polycone(shape, point, dir)
    elseif shape isa CutTube
        distanceToIn_cuttube(shape, point, dir)
    elseif shape isa BooleanUnion
        distanceToIn_booleanunion(shape, point, dir)
    elseif shape isa BooleanSubtraction
        distanceToIn_booleansubtraction(shape, point, dir)
    elseif shape isa BooleanIntersection
        distanceToIn_booleanintersection(shape, point, dir)
    elseif shape isa PlacedVolume
        xf = shape.transformation
        distanceToIn_placedvolume(shape.volume.shape, xf*point, xf*dir)
    end
end
function distanceToOut(shape, point, dir)
    if shape isa Trap
        distanceToOut_trap(shape, point, dir)
    elseif shape isa Trd
        distanceToOut_trd(shape, point, dir)
    elseif shape isa Cone
        distanceToOut_cone(shape, point, dir)
    elseif shape isa Box
        distanceToOut_box(shape, point, dir)
    elseif shape isa Tube
        distanceToOut_tube(shape, point, dir)
    elseif shape isa Polycone
        distanceToOut_polycone(shape, point, dir)
    elseif shape isa CutTube
        distanceToOut_cuttube(shape, point, dir)
    elseif shape isa BooleanUnion
        distanceToOut_booleanunion(shape, point, dir)
    elseif shape isa BooleanSubtraction
        distanceToOut_booleansubtraction(shape, point, dir)
    elseif shape isa BooleanIntersection
        distanceToOut_booleanintersection(shape, point, dir)
    end
end
