function distanceToIn(shape, point, dir)
    if shape isa Trap
        distanceToIn_trap(shape::Trap, point, dir)
    elseif shape isa Trd
        distanceToIn_trd(shape::Trd, point, dir)
    elseif shape isa Cone
        distanceToIn_cone(shape::Cone, point, dir)
    elseif shape isa Box
        distanceToIn_box(shape::Box, point, dir)
    elseif shape isa Tube
        distanceToIn_tube(shape::Tube, point, dir)
    elseif shape isa Polycone
        distanceToIn_polycone(shape::Polycone, point, dir)
    elseif shape isa CutTube
        distanceToIn_cuttube(shape::CutTube, point, dir)
    elseif shape isa BooleanUnion
        distanceToIn_booleanunion(shape::BooleanUnion, point, dir)
    elseif shape isa BooleanSubtraction
        distanceToIn_booleansubtraction(shape::BooleanSubtraction, point, dir)
    elseif shape isa BooleanIntersection
        distanceToIn_booleanintersection(shape::BooleanIntersection, point, dir)
    end
end
function distanceToOut(shape, point, dir)
    if shape isa Trap
        distanceToOut_trap(shape::Trap, point, dir)
    elseif shape isa Trd
        distanceToOut_trd(shape::Trd, point, dir)
    elseif shape isa Cone
        distanceToOut_cone(shape::Cone, point, dir)
    elseif shape isa Box
        distanceToOut_box(shape::Box, point, dir)
    elseif shape isa Tube
        distanceToOut_tube(shape::Tube, point, dir)
    elseif shape isa Polycone
        distanceToOut_polycone(shape::Polycone, point, dir)
    elseif shape isa CutTube
        distanceToOut_cuttube(shape::CutTube, point, dir)
    elseif shape isa BooleanUnion
        distanceToOut_booleanunion(shape::BooleanUnion, point, dir)
    elseif shape isa BooleanSubtraction
        distanceToOut_booleansubtraction(shape::BooleanSubtraction, point, dir)
    elseif shape isa BooleanIntersection
        distanceToOut_booleanintersection(shape::BooleanIntersection, point, dir)
    end
end
