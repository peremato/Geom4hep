using Test
using Geom4hep, LinearAlgebra

@testset "Geom4hep tests" verbose = true begin 
    include("testTransformation3D.jl")
    include("testBox.jl")
    include("testTrd.jl")
    include("testTube.jl")
    include("testCone.jl")
    include("testPolycone.jl")
    include("testPolyhedra.jl")
    include("testTorus.jl")
    include("testTrap.jl")
    include("testCutTube.jl")
    include("testSphere.jl")
    include("testBoolean.jl")
    include("testMultiUnion.jl")
    include("testVolume.jl")
    include("testNormals.jl")
    include("testNavigation.jl")
    include("testModel.jl")
    include("testGDML.jl")
    include("testBenchmark.jl")
end
