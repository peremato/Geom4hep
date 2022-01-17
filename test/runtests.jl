using Test
using Geom4hep, LinearAlgebra

@testset "Geom4hep tests" verbose = true begin 
    include("testTransformation3D.jl")
    include("testBox.jl")
    include("testTrd.jl")
    include("testTube.jl")
    include("testCone.jl")
end

