@testset "Box" begin
    pzero = Point3(0., 0., 0.)
    ponxside, ponyside, ponzside = Point3(20., 0., 0.), Point3(0., 30., 0.), Point3(0., 0., 40.)
    ponmxside, ponmyside, ponmzside = Point3(-20., 0., 0.), Point3(0., -30., 0.), Point3(0., 0., -40.)
    ponzsidey, ponmzsidey = Point3(0., 25., 40.), Point3(0., 25., -40.)
    
    pbigx = Point3(100., .0, 0.)
    pbigy = Point3(0., 100., 0.)
    pbigz = Point3(0., 0., 100.)
    pbigmx = Point3(-100., 0., 0.)
    pbigmy = Point3(0., -100., 0.)
    pbigmz = Point3(0., 0., -100.)
    
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
    @test minExtent ≈ Point3(-20, -30, -40)
    @test maxExtent ≈ Point3(20, 30, 40)
    minExtent, maxExtent = extent(b2)
    @test minExtent ≈ Point3(-10, -10, -10)
    @test maxExtent ≈ Point3(10, 10, 10)

    # Check Surface Normal
    nvect = normal(b1, ponxside)
    @test normal(b1, ponxside) ≈ Vector3(1,0,0)
    @test normal(b1, ponmxside) ≈ Vector3(-1,0,0)
    @test normal(b1, ponyside) ≈ Vector3(0,1,0)
    @test normal(b1, ponmyside) ≈ Vector3(0,-1,0)
    @test normal(b1, ponzside) ≈ Vector3(0,0,1)
    @test normal(b1, ponmzside) ≈ Vector3(0,0,-1)
    @test normal(b1, ponzsidey) ≈ Vector3(0,0,1)
    @test normal(b1, ponmzsidey) ≈ Vector3(0,0,-1)

    # Normals on Edges
    edgeXY = Point3(20.0, 30., 0.0)
    edgemXmY = Point3(-20.0, -30., 0.0)
    edgeXmY = Point3(20.0, -30., 0.0)
    edgemXY = Point3(-20.0, 30., 0.0)
    edgeXZ = Point3(20.0, 0.0, 40.0)
    edgemXmZ = Point3(-20.0, 0.0, -40.0)
    edgeXmZ = Point3(20.0, 0.0, -40.0)
    edgemXZ = Point3(-20.0, 0.0, 40.0)
    edgeYZ = Point3(0.0, 30.0, 40.0)
    edgemYmZ = Point3(0.0, -30.0, -40.0)
    edgeYmZ = Point3(0.0, 30.0, -40.0)
    edgemYZ = Point3(0.0, -30.0, 40.0)

    @test normal(b1, edgeXY) ≈ Vector3(1/√2, 1/√2, 0.0)
    @test normal(b1, edgemXmY) ≈ Vector3(-1/√2, -1/√2, 0.0)
    @test normal(b1, edgeXmY) ≈ Vector3(1/√2, -1/√2, 0.0)
    @test normal(b1, edgemXY) ≈ Vector3(-1/√2, 1/√2, 0.0)

    @test normal(b1, edgeXZ) ≈ Vector3(1/√2, 0.0, 1/√2)
    @test normal(b1, edgemXmZ) ≈ Vector3(-1/√2, 0.0, -1/√2)
    @test normal(b1, edgeXmZ) ≈ Vector3(1/√2, 0.0, -1/√2)
    @test normal(b1, edgemXZ) ≈ Vector3(-1/√2, 0.0, 1/√2)

    @test normal(b1, edgeYZ) ≈ Vector3(0.0, 1/√2, 1/√2)
    @test normal(b1, edgemYmZ) ≈ Vector3(0.0, -1/√2, -1/√2)
    @test normal(b1, edgeYmZ) ≈ Vector3(0.0, 1/√2, -1/√2)
    @test normal(b1, edgemYZ) ≈ Vector3(0.0, -1/√2, 1/√2)


    # Normals on corners
    cornerXYZ = Point3(20.0, 30., 40.0)
    cornermXYZ = Point3(-20.0, 30., 40.0)
    cornerXmYZ = Point3(20.0, -30., 40.0)
    cornermXmYZ = Point3(-20.0, -30., 40.0)
    cornerXYmZ = Point3(20.0, 30., -40.0)
    cornermXYmZ = Point3(-20.0, 30., -40.0)
    cornerXmYmZ = Point3(20.0, -30., -40.0)
    cornermXmYmZ = Point3(-20.0, -30., -40.0)

    @test normal(b1, cornerXYZ) ≈ Vector3(1/√3, 1/√3, 1/√3)
    @test normal(b1, cornermXYZ) ≈ Vector3(-1/√3, 1/√3, 1/√3)
    @test normal(b1, cornerXmYZ) ≈ Vector3(1/√3, -1/√3, 1/√3)
    @test normal(b1, cornermXmYZ) ≈ Vector3(-1/√3, -1/√3, 1/√3)
    @test normal(b1, cornerXYmZ) ≈ Vector3(1/√3, 1/√3, -1/√3)
    @test normal(b1, cornermXYmZ) ≈ Vector3(-1/√3, 1/√3, -1/√3)
    @test normal(b1, cornerXmYmZ) ≈ Vector3(1/√3, -1/√3, -1/√3)
    @test normal(b1, cornermXmYmZ) ≈ Vector3(-1/√3, -1/√3, -1/√3)

    vx = Vector3(1., 0., 0.)
    vy = Vector3(0., 1., 0.)
    vz = Vector3(0., 0., 1.)
    vmx = Vector3(-1., 0., 0.)
    vmy = Vector3(0., -1., 0.)
    vmz = Vector3(0., 0., -1.)
    vxy = Vector3(1/√2, 1/√2, 0.)
    vmxy = Vector3(-1/√2, 1/√2, 0.)
    vmxmy = Vector3(-1/√2, -1/√2, 0.)
    vxmy = Vector3(1/√2, -1/√2, 0.)
    vxmz = Vector3(1/√2, 0., -1/√2)

    # distanceToOut(P,V) with asserts for normal and convex
    dist = distanceToOut(b1, pzero, vx)
    norm = normal(b1, pzero + dist * vx) 
    @test dist ≈ 20 && norm ≈ vx

    dist = distanceToOut(b1, pzero, vmx)
    norm = normal(b1, pzero + dist * vmx) 
    @test dist ≈ 20 && norm ≈ vmx 

    dist = distanceToOut(b1, pzero, vy)
    norm = normal(b1, pzero + dist * vy) 
    @test dist ≈ 30 && norm ≈ vy 

    dist = distanceToOut(b1, pzero, vmy)
    norm = normal(b1, pzero + dist * vmy) 
    @test dist ≈ 30 && norm ≈ vmy 

    dist = distanceToOut(b1, pzero, vz)
    norm = normal(b1, pzero + dist * vz) 
    @test dist ≈ 40 && norm ≈ vz 

    dist = distanceToOut(b1, pzero, vmz)
    norm = normal(b1, pzero + dist * vmz) 
    @test dist ≈ 40 && norm ≈ vmz 

    dist = distanceToOut(b1, pzero, vxy)
    norm = normal(b1, pzero + dist * vxy) 
    @test dist ≈ √800


    #testing a few special directions ( with special zero values )
    # important since we sometimes operate with sign bits etc.
    @test distanceToOut(b1, pzero, Vector3(-0., -0., -1.)) ≈ 40.
    @test distanceToOut(b1, pzero, Vector3( 0., -0., -1.)) ≈ 40.
    @test distanceToOut(b1, pzero, Vector3( 0.,  0., -1.)) ≈ 40.
    @test distanceToOut(b1, pzero, Vector3( 0.,  1.,  0.)) ≈ 30.
    @test distanceToOut(b1, pzero, Vector3( 0.,  1., -0.)) ≈ 30.
    @test distanceToOut(b1, pzero, Vector3(-0.,  1.,  0.)) ≈ 30.
    @test distanceToOut(b1, pzero, Vector3( 1., -0.,  0.)) ≈ 20.
    @test distanceToOut(b1, pzero, Vector3( 1.,  0., -0.)) ≈ 20.
    @test distanceToOut(b1, pzero, Vector3( 1., -0., -0.)) ≈ 20.
    
    dist = distanceToOut(b1, ponxside, vx)
    norm = normal(b1, ponxside + dist * vx)
    @test dist ≈ 0. && norm ≈ vx
    
    dist = distanceToOut(b1, ponmxside, vmx)
    norm = normal(b1, ponmxside + dist * vmx)
    @test dist ≈ 0. && norm ≈ vmx
    
    dist = distanceToOut(b1, ponmxside, vmx)
    norm = normal(b1, ponmxside + dist * vmx)
    @test dist ≈ 0. && norm ≈ vmx
    
    dist = distanceToOut(b1, ponyside, vy)
    norm = normal(b1, ponyside + dist * vy)
    @test dist ≈ 0. && norm ≈ vy
    
    dist = distanceToOut(b1, ponmyside, vmy)
    norm = normal(b1, ponmyside + dist * vmy)
    @test dist ≈ 0. && norm ≈ vmy
    
    dist = distanceToOut(b1, ponzside, vz)
    norm = normal(b1, ponzside + dist * vz)
    @test dist ≈ 0. && norm ≈ vz
    
    dist = distanceToOut(b1, ponmzside, vmz)
    norm = normal(b1, ponmzside + dist * vmz)
    @test dist ≈ 0. && norm ≈ vmz

    # Check Inside
    @test inside(b1, pzero) == kInside
    @test inside(b1, pbigz) == kOutside
    @test inside(b1, ponxside) == kSurface
    @test inside(b1, ponyside) == kSurface
    @test inside(b1, ponzside) == kSurface

    @test inside(b2, pzero) == kInside
    @test inside(b2, pbigz) == kOutside
    @test inside(b2, ponxside) == kOutside
    @test inside(b2, ponyside) == kOutside
    @test inside(b2, ponzside) == kOutside
    @test inside(b2, Point3(10., 0., 0.)) == kSurface
    @test inside(b2, Point3(0., 10., 0.)) == kSurface
    @test inside(b2, Point3(0., 0., 10.)) == kSurface
    @test inside(b2, Point3(20., 20., 10.)) == kOutside
    @test inside(b2, Point3(100., 10., 30.)) == kOutside
    @test inside(b2, Point3(10., 20., 20.)) == kOutside

    # SafetyToOut(P)
    @test safetyToOut(b1, pzero) ≈ 20
    @test safetyToOut(b1, Point3(1.,0.,0.)) ≈ 19
    @test safetyToOut(b1, Point3(0.,1.,0.)) ≈ 20
    @test safetyToOut(b1, Point3(0.,0.,1.)) ≈ 20

    # SafetyToIn(P)
    @test safetyToIn(b1, pbigx) ≈ 80
    @test safetyToIn(b1, pbigmx) ≈ 80
    @test safetyToIn(b1, pbigy) ≈ 70
    @test safetyToIn(b1, pbigmy) ≈ 70
    @test safetyToIn(b1, pbigz) ≈ 60
    @test safetyToIn(b1, pbigmz) ≈ 60

    # DistanceToIn(P,V)
    @test distanceToIn(b1, pbigx, vmx) ≈ 80
    @test distanceToIn(b1, pbigmx, vx) ≈ 80
    @test distanceToIn(b1, pbigy, vmy) ≈ 70
    @test distanceToIn(b1, pbigmy, vy) ≈ 70
    @test distanceToIn(b1, pbigz, vmz) ≈ 60
    @test distanceToIn(b1, pbigmz, vz) ≈ 60

    # testing a few special directions ( with special zero values )
    # important since we sometimes operate with sign bits etc.
    @test distanceToIn(b1, Point3(0., 0., -50.), Vector3(-0., 0., 1.)) ≈ 10
    @test distanceToIn(b1, Point3(0., 0., -50.), Vector3( 0.,-0., 1.)) ≈ 10
    @test distanceToIn(b1, Point3(0., 0.,  50.), Vector3( 0., 0.,-1.)) ≈ 10
    @test distanceToIn(b1, Point3(0., -40., 0.), Vector3( 0., 1., 0.)) ≈ 10
    @test distanceToIn(b1, Point3(0., -40., 0.), Vector3( 0., 1.,-0.)) ≈ 10
    @test distanceToIn(b1, Point3(0., -40., 0.), Vector3(-0., 1., 0.)) ≈ 10
    @test distanceToIn(b1, Point3(-30., 0., 0.), Vector3( 1.,-0., 0.)) ≈ 10
    @test distanceToIn(b1, Point3( 30,  0., 0.), Vector3(-1., 0.,-0.)) ≈ 10
    @test distanceToIn(b1, Point3(-30., 0., 0.), Vector3( 1.,-0.,-0.)) ≈ 10

    @test distanceToIn(b1, pbigx, vxy) ≈ Inf
    @test distanceToIn(b1, pbigmx, vxy) ≈ Inf

    pJohnXZ = Point3(9., 0., 12.)
    @test distanceToIn(b2, pJohnXZ, vxmz) ≈ Inf
    pJohnXY = Point3(12., 9., 0.)
    @test distanceToIn(b2, pJohnXY, vxmy) ≈ Inf
    @test distanceToIn(b2, pJohnXY, vmx) ≈ 2
    pMyXY = Point3(32., -11., 0.)
    @test distanceToIn(b2, pMyXY, vmxy) ≈ Inf
    @test distanceToIn(b1, Point3(-25., -35., 0.), vx) ≈ Inf
    @test distanceToIn(b1, Point3(-25., -35., 0.), vy) ≈ Inf
    @test distanceToIn(b2, pJohnXY, vmx) ≈ 2

    tempDir = normalize(Vector3(-0.76165597579890043, 0.64364445891356026, -0.074515708658524193))
    d = distanceToIn(box3, Point3(0.15000000185, -22.048743592955137, 2.4268539333219472), tempDir)
    @test isapprox(d, 0, atol=1e-6)

    # testing tolerance of DistanceToIn
    b4 = Box(5., 5., 5.)
    # a point very slightly inside should return 0
    tempDir = normalize(Vector3(0.76315134679548990437, 0.53698876104646497964, -0.35950395323836459305))
    @test distanceToIn(b4, Point3(-3.0087437277453119577, -5.0+kTolerance()/2, 4.8935648380409944025), tempDir) <= 0.0


    # a point on the surface pointing outside must return infinity length
    @test distanceToIn(b2,Point3(10., 0., 0.), Vector3(1., 0., 0.)) == Inf
    @test distanceToIn(b2,Point3(10. - 0.9*eps(), 0., 0.), Vector3(1.,0.,0.)) >= Inf

end

@testset "Box_VG_431" begin
    vx = Vector3(1., 0., 0.)
    vy = Vector3(0., 1., 0.)
    vz = Vector3(0., 0., 1.)
  
    bx = Box(1., 200., 200.)
    by = Box(200., 1., 200.)
    bz = Box(200., 200., 1.)
    # slightly outside of +x face and entering and exiting
    testp = Point3(1.,0.,0.) + 6e-15 * vx
    testv = normalize(Vector3(0.,1.,1.) - 4e-6 * vx)
    dist = distanceToIn(bx, testp, testv)
    @test isapprox(dist, 0, atol=1e-6)

    dist = distanceToOut(bx, testp, testv)
    norm = normal(bx, testp + dist*testv)
    @test dist ≈ 200*√2
    @test norm ≈ normalize(vy+vz)

    dist = distanceToOut(bx, testp, vx)
    norm = normal(bx, testp + dist * vx)
    @test isapprox(dist, 0., atol=1e-6)
    @test norm ≈ vx

    # slightly outside of -x face and entering and exiting
    testp = Point3(-1., 0., 0.) - 6e-15 * vx
    testv = normalize(Vector3(0., 1., 1.) + 4e-6 * vx)
    dist = distanceToIn(bx, testp, testv)
    @test isapprox(dist, 0, atol=1e-6)

    dist = distanceToOut(bx, testp, testv)
    norm = normal(bx, testp + dist*testv)
    @test dist ≈ 200*√2
    @test norm ≈ normalize(vy+vz)

    dist = distanceToOut(bx, testp, -vx)
    norm = normal(bx, testp - dist * vx)
    @test isapprox(dist, 0., atol=1e-6)
    @test norm ≈ -vx

    # slightly outside of +y face and entering and exiting
    testp = Point3(0.,1.,0.) + 6e-15 * vy
    testv = normalize(Vector3(1.,0.,1.) - 4e-6 * vy)
    dist = distanceToIn(by, testp, testv)
    @test isapprox(dist, 0, atol=1e-6)

    dist = distanceToOut(by, testp, testv)
    norm = normal(by, testp + dist*testv)
    @test dist ≈ 200*√2
    @test norm ≈ normalize(vx+vz)

    dist = distanceToOut(by, testp, vy)
    norm = normal(by, testp + dist * vy)
    @test isapprox(dist, 0., atol=1e-6)
    @test norm ≈ vy

    # slightly outside of -y face and entering and exiting
    testp = Point3(0.,-1.,0.) - 6e-15 * vy
    testv = normalize(Vector3(1.,0.,1.) + 4e-6 * vy)
    dist = distanceToIn(by, testp, testv)
    @test isapprox(dist, 0, atol=1e-6)

    dist = distanceToOut(by, testp, testv)
    norm = normal(by, testp + dist*testv)
    @test dist ≈ 200*√2
    @test norm ≈ normalize(vx+vz)

    dist = distanceToOut(by, testp, -vy)
    norm = normal(by, testp - dist * vy)
    @test isapprox(dist, 0., atol=1e-6)
    @test norm ≈ -vy

    # slightly outside of +z face and entering and exiting
    testp = Point3(0.,0.,1.) + 6e-15 * vz
    testv = normalize(Vector3(1.,1., 0.) - 4e-6 * vz)
    dist = distanceToIn(bz, testp, testv)
    @test isapprox(dist, 0, atol=1e-6)

    dist = distanceToOut(bz, testp, testv)
    norm = normal(bz, testp + dist*testv)
    @test dist ≈ 200*√2
    @test norm ≈ normalize(vx+vy)

    dist = distanceToOut(bz, testp, vz)
    norm = normal(bz, testp + dist * vz)
    @test isapprox(dist, 0., atol=1e-6)
    @test norm ≈ vz

    # slightly outside of -z face and entering and exiting
    testp = Point3(0.,0.,-1.) - 6e-15 * vz
    testv = normalize(Vector3(1.,1.,0.) + 4e-6 * vz)
    dist = distanceToIn(bz, testp, testv)
    @test isapprox(dist, 0, atol=1e-6)

    dist = distanceToOut(bz, testp, testv)
    norm = normal(bz, testp + dist*testv)
    @test dist ≈ 200*√2
    @test norm ≈ normalize(vx+vy)

    dist = distanceToOut(bz, testp, -vz)
    norm = normal(bz, testp - dist * vz)
    @test isapprox(dist, 0., atol=1e-6)
    @test norm ≈ -vz

    # slightly inside of +x face
    testp = Point3(1.,0.,0.) - 6e-15 * vx
    testv = normalize(Vector3(0.,1.,1.) + 4e-6 * vx)
    dist = distanceToIn(bx, testp, testv)
    @test dist ≈  Inf
    dist = distanceToOut(bx, testp, testv)
    @test isapprox(dist, 0., atol=1e-6)

    # slightly inside of -x face
    testp = Point3(-1.,0.,0.) + 6e-15 * vx
    testv = normalize(Vector3(0.,1.,1.) - 4e-6 * vx)
    dist = distanceToIn(bx, testp, testv)
    @test dist ≈  Inf
    dist = distanceToOut(bx, testp, testv)
    @test isapprox(dist, 0., atol=1e-6)

    # slightly outside of +x face and exiting
    cube = Box(10., 10., 10.)
    testp = Point3(10., 0., 0.) + 4.e-10 * vx
    dist = distanceToOut(cube, testp, vx)
    @test isapprox(dist, 0., atol=1e-6)
    norm = normal(cube, testp + dist * vz)
    @test norm ≈ vx
end
