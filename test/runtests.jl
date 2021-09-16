using Test
using Geom4hep

@testset "Box" begin
    pzero = Point(0., 0., 0.)
    ponxside, ponyside, ponzside = Point(20., 0., 0.), Point(0., 30., 0.), Point(0., 0., 40.)
    ponmxside, ponmyside, ponmzside = Point(-20., 0., 0.), Point(0., -30., 0.), Point(0., 0., -40.)
    ponzsidey, ponmzsidey = Point(0., 25., 40.), Point(0., 25., -40.)
    
    #  pbigx(100, 0, 0), pbigy(0, 100, 0), pbigz(0, 0, 100);
    #  pbigmx(-100, 0, 0), pbigmy(0, -100, 0), pbigmz(0, 0, -100);
    #  vx(1, 0, 0), vy(0, 1, 0), vz(0, 0, 1);
    #  vmx(-1, 0, 0), vmy(0, -1, 0), vmz(0, 0, -1);
    #  vxy(1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
    #  vmxy(-1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
    #  vmxmy(-1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
    #  vxmy(1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
    #  vxmz(1 / std::sqrt(2.0), 0, -1 / std::sqrt(2.0));

    b1 = Box(20., 30., 40.)
    b2 = Box(10., 10., 10.)
    box3 = Box(0.14999999999999999, 24.707000000000001, 22.699999999999999)

    # Volume and Surfaces
    @test capacity(b2) == 8000
    @test capacity(b1) == 192000
    @test surface(b1) == 20800
    @test surface(b2) == 6 * 20 * 20

    # CalculateExtent
    minExtent, maxExtent = extent(b1)
    @test minExtent ≈ Point(-20, -30, -40)
    @test maxExtent ≈ Point(20, 30, 40)
    minExtent, maxExtent = extent(b2)
    @test minExtent ≈ Point(-10, -10, -10)
    @test maxExtent ≈ Point(10, 10, 10)

    # Check Surface Normal
    nvect = normal(b1, ponxside)
    @test normal(b1, ponxside) ≈ Vec(1,0,0)
    @test normal(b1, ponmxside) ≈ Vec(-1,0,0)
    @test normal(b1, ponyside) ≈ Vec(0,1,0)
    @test normal(b1, ponmyside) ≈ Vec(0,-1,0)
    @test normal(b1, ponzside) ≈ Vec(0,0,1)
    @test normal(b1, ponmzside) ≈ Vec(0,0,-1)
    @test normal(b1, ponzsidey) ≈ Vec(0,0,1)
    @test normal(b1, ponmzsidey) ≈ Vec(0,0,-1)

    # Normals on Edges
    edgeXY = Point(20.0, 30., 0.0)
    edgemXmY = Point(-20.0, -30., 0.0)
    edgeXmY = Point(20.0, -30., 0.0)
    edgemXY = Point(-20.0, 30., 0.0)
    edgeXZ = Point(20.0, 0.0, 40.0)
    edgemXmZ = Point(-20.0, 0.0, -40.0)
    edgeXmZ = Point(20.0, 0.0, -40.0)
    edgemXZ = Point(-20.0, 0.0, 40.0)
    edgeYZ = Point(0.0, 30.0, 40.0)
    edgemYmZ = Point(0.0, -30.0, -40.0)
    edgeYmZ = Point(0.0, 30.0, -40.0)
    edgemYZ = Point(0.0, -30.0, 40.0)

    @test normal(b1, edgeXY) ≈ Vec(1/√2, 1/√2, 0.0)
    @test normal(b1, edgemXmY) ≈ Vec(-1/√2, -1/√2, 0.0)
    @test normal(b1, edgeXmY) ≈ Vec(1/√2, -1/√2, 0.0)
    @test normal(b1, edgemXY) ≈ Vec(-1/√2, 1/√2, 0.0)

    @test normal(b1, edgeXZ) ≈ Vec(1/√2, 0.0, 1/√2)
    @test normal(b1, edgemXmZ) ≈ Vec(-1/√2, 0.0, -1/√2)
    @test normal(b1, edgeXmZ) ≈ Vec(1/√2, 0.0, -1/√2)
    @test normal(b1, edgemXZ) ≈ Vec(-1/√2, 0.0, 1/√2)

    @test normal(b1, edgeYZ) ≈ Vec(0.0, 1/√2, 1/√2)
    @test normal(b1, edgemYmZ) ≈ Vec(0.0, -1/√2, -1/√2)
    @test normal(b1, edgeYmZ) ≈ Vec(0.0, 1/√2, -1/√2)
    @test normal(b1, edgemYZ) ≈ Vec(0.0, -1/√2, 1/√2)


    # Normals on corners
    cornerXYZ = Point(20.0, 30., 40.0)
    cornermXYZ = Point(-20.0, 30., 40.0)
    cornerXmYZ = Point(20.0, -30., 40.0)
    cornermXmYZ = Point(-20.0, -30., 40.0)
    cornerXYmZ = Point(20.0, 30., -40.0)
    cornermXYmZ = Point(-20.0, 30., -40.0)
    cornerXmYmZ = Point(20.0, -30., -40.0)
    cornermXmYmZ = Point(-20.0, -30., -40.0)

    @test normal(b1, cornerXYZ) ≈ Vec(1/√3, 1/√3, 1/√3)
    @test normal(b1, cornermXYZ) ≈ Vec(-1/√3, 1/√3, 1/√3)
    @test normal(b1, cornerXmYZ) ≈ Vec(1/√3, -1/√3, 1/√3)
    @test normal(b1, cornermXmYZ) ≈ Vec(-1/√3, -1/√3, 1/√3)
    @test normal(b1, cornerXYmZ) ≈ Vec(1/√3, 1/√3, -1/√3)
    @test normal(b1, cornermXYmZ) ≈ Vec(-1/√3, 1/√3, -1/√3)
    @test normal(b1, cornerXmYmZ) ≈ Vec(1/√3, -1/√3, -1/√3)
    @test normal(b1, cornermXmYmZ) ≈ Vec(-1/√3, -1/√3, -1/√3)

end
