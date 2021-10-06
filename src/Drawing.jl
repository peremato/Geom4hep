using GLMakie
using Colors

colors = colormap("Grays", 3)
#---Draw a Volume---------------------------------------------------------------
function draw(s::Scene, vol::Volume, t::Transformation3D, level::Int64)
    m = toMesh(vol.shape)
    if isone(t)
        mesh!(s, m, color=colors[level], transparency=true, ambient=0.7, visible= vol.label == "world" ? false : true)
    else
        points = GeometryBasics.coordinates(m)
        faces  = GeometryBasics.faces(m)
        map!(c -> c * t, points, points)
        mesh!(s, points, faces, color=colors[level], transparency=true, ambient=0.7)
    end
    for daughter in vol.daughters
        draw(s, daughter.volume, daughter.transformation * t, level+1)
    end
end

function draw(s::Scene, vol::Volume)
    draw(s, vol, one(Transformation3D{Float64}), 1)
    display(s)
end

