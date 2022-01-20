using Test
using Geom4hep, LinearAlgebra

include("../src/Units.jl")

@testset "Cone{$T}" for T in (Float64, Float32)
    # Cones
    c1  = Cone{T}(50, 100, 50, 100, 50, 0, 2π)
    c2  = Cone{T}(50, 100, 50, 200, 50, -1, 2π)
    c3  = Cone{T}(50, 100, 50, 100, 50, -π/6, π/3)
    c4  = Cone{T}(50, 100, 50, 200, 50, -π/6, π/3)
    c5  = Cone{T}(25, 50, 75, 150, 50, 0, 3π/2)
    c6  = Cone{T}(0, 150, 0, 150, 50, 0, 2π)
    c7  = Cone{T}(95, 100, 95, 100, 50, 0, 2π)
    c8a = Cone{T}(0, 100, 0, 150, 50, 0, 2π)
    c8b = Cone{T}(50, 100, 100, 150, 50, 0, 2π)
    c8c = Cone{T}(100, 150, 50, 100, 50, 0, 2π)
    ctest10 = Cone{T}(20, 60, 80, 140, 100, 10deg, 300deg)
    c9  = Cone{T}(50, 60, 0, 10, 50, 0, 2π)
    cms = Cone{T}(0.0, 70.0, 0.0, 157.8, 2949.0, 0.0, 2π)
    cms2= Cone{T}(401.0, 1450.0, 1020.0, 1450.0, 175.0, 0.0, 2π)

    # points
    pzero = Point3{T}(0, 0, 0)
    pplx = Point3{T}(120, 0, 0)
    pply = Point3{T}(0, 120, 0)
    pplz = Point3{T}(0, 0, 120)
    pmix = Point3{T}(-120, 0, 0)
    pmiy = Point3{T}(0, -120, 0)
    pmiz = Point3{T}(0, 0, -120)
    ponmiz = Point3{T}(0, 75, -50)
    ponplz = Point3{T}(0, 75, 50)
    ponr1 = Point3{T}(50/√2, 50/√2, 0)
    ponr2 = Point3{T}(100/√2, 100/√2, 0)
    ponphi1 = Point3{T}(60*cos(π/6), -60*sin(π/6), 0)
    ponphi2 = Point3{T}(60*cos(π/6),  60*sin(π/6), 0)
    ponr2b = Point3{T}(150, 0, 0)
    pnearplz = Point3{T}(45, 45, 45)
    pnearmiz = Point3{T}(45, 45, -45)
    pydx = Point3{T}(60, 150, 0)
    pbigx = Point3{T}(500, 0, 0)
    proot1 = Point3{T}(0, 125, -1000)
    proot2 = Point3{T}(0, 75, -1000)
    pparr1 = Point3{T}(0, 25, -150) # Test case parallel to both rs of c8
    pparr2 = Point3{T}(0, 75, -50)
    pparr3 = Point3{T}(0, 125, 50)
    vparr  = Vector3{T}(0, 1/√5, 2/√5)
    vnphi1 = Vector3{T}(-sin(π/6),-cos(π/6), 0)
    vnphi2 = Vector3{T}(-sin(π/6), cos(π/6), 0)
    pct10   =   Point3{T}(60, 0, 0)
    pct10mx =   Point3{T}(-50, 0, 0)
    pct10phi1 = Point3{T}(60*cos(10deg), 60*sin(10deg), 0)
    pct10phi2 = Point3{T}(60*cos(50deg),-60*sin(50deg), 0)
    pct10e1 = Point3{T}(-691 - 500, 174, 404)
    pct10e2 = Point3{T}(400 - 500, 20.9, 5.89)
    pct10e3 = Point3{T}(456 - 500, 13, -14.7)
    pct10e4 = Point3{T}(537 - 500, 1.67, -44.1)
    # point P is outside
    pct10e5 = Point3{T}(537, 1.67, -44.1)
    pct10e6 = Point3{T}(1e+03, -63.5, -213)
    pt10s1  = Point3{T}(6.454731216775542, -90.42080754048007, 100.)
    pt10s2  = Point3{T}(22.65282328600368, -69.34877585931267, 76.51600623610082)
    pt10s3  = Point3{T}(51.28206938732319, -32.10510677306267, 35.00932544708616)
    vt10d   = Vector3{T}(0.4567090876640433, 0.5941309830320264, -0.6621368319663807)


    vx    = Vector3{T}( 1, 0, 0)
    vy    = Vector3{T}( 0, 1, 0)
    vz    = Vector3{T}( 0, 0, 1)
    vmx   = Vector3{T}(-1, 0, 0)    
    vmy   = Vector3{T}( 0,-1, 0)
    vmz   = Vector3{T}( 0, 0,-1)
    vxy   = Vector3{T}( 1/√2, 1/√2,0)
    vmxy  = Vector3{T}(-1/√2, 1/√2,0)
    vxmy  = Vector3{T}( 1/√2,-1/√2,0)
    vmxmy = Vector3{T}(-1/√2,-1/√2,0)
    vx2mz = Vector3{T}( 1/√5, 0, -2/√5)
    vxmz  = Vector3{T}( 1/√2, 0, -1/√2)

    d1 = normalize(pct10e2 - pct10e1)

    @test capacity(c1) ≈ T(2π) * 50 * (100 * 100 - 50 * 50)
    @test capacity(c6) ≈ T(2π) * 50 * (150 * 150)

    #@test surface(c1) ≈ T(2π) * (50 * 2 * 50 + 100 * 2 * 50 + 100 * 100 - 50 * 50)

    @test inside(ctest10, pct10e1) == kOutside
    @test inside(ctest10, pct10e2) == kInside
    @test inside(ctest10, pct10e3) == kInside
    @test inside(ctest10, pct10e4) == kOutside
    @test inside(ctest10, pct10e5) == kOutside
    @test inside(ctest10, pct10e6) == kOutside
    @test inside(ctest10, pct10mx) == kSurface
    @test inside(ctest10, pt10s1) == kSurface
    @test inside(ctest10, pt10s2) == kSurface
    @test inside(ctest10, pt10s3) == kOutside

    @test inside(c1, pzero) == kOutside
    @test inside(c6, pzero) == kInside
    @test inside(c1, pplx) == kOutside
    @test inside(c2, pplx) == kInside
    @test inside(c3, pplx) == kOutside
    @test inside(c4, pplx) == kInside
    @test inside(c1, ponmiz) == kSurface
    @test inside(c1, ponplz) == kSurface
    @test inside(c1, ponr1) == kSurface
    @test inside(c1, ponmiz) == kSurface
    @test inside(c1, ponr2) == kSurface
    @test inside(c3, ponphi1) == kSurface
    @test inside(c3, ponphi2) == kSurface
    @test inside(c5, Point3{T}(70,1,0)) == kInside
    @test inside(c5, Point3{T}(50,-50,0)) == kOutside
    @test inside(c5, Point3{T}(70,0,0)) == kSurface
    @test inside(c5, Point3{T}(100,0,0)) == kSurface
    @test inside(c3, Point3{T}(100,0,0)) == kSurface
    @test inside(c5, Point3{T}(100,0,50)) == kSurface
    @test inside(c3, Point3{T}(100,0,50)) == kSurface

    # safety
    @test safetyToOut(c4, ponphi1) ≈ 0.
    @test safetyToOut(c1, ponphi1) ≈ 10.
    @test safetyToOut(c1, pnearplz) ≈ 5.
    @test safetyToOut(c1, pnearmiz) ≈ 5.
    @test isapprox(safetyToOut(c1, ponr1), 0., atol=100*eps(T))
    @test isapprox(safetyToOut(c1, ponr2), 0., atol=100*eps(T))
    @test safetyToOut(c6, pzero) ≈ 50.
    @test isapprox(safetyToOut(c5, Point3{T}(0,-70, 0)), 0., atol=100*eps(T))

    # distToOut
    @test distanceToOut(c4, pplx, vx) ≈ 30.
    @test distanceToOut(c2, pplx, vx) ≈ 30.
    @test distanceToOut(c4, pplx, vmx) ≈ 70.
    @test distanceToOut(c2, pplx, vmx) ≈ 70.
    @test distanceToOut(c3, ponphi1, vmy) ≈ 0.
    @test distanceToOut(c3, ponphi1, vy) ≈ 120 * sin(π/6)
    @test distanceToOut(c3, ponphi2, vy) ≈ 0.
    @test distanceToOut(c3, ponphi2, vmy) ≈ 120 * sin(π/6)
    @test distanceToOut(c6, ponplz, vmz) ≈ 100.
    @test distanceToOut(c6, ponplz, vz) ≈ 0.
    @test distanceToOut(c6, ponmiz, vz) ≈ 100.
    @test distanceToOut(c6, ponmiz, vmz) ≈ 0.
    @test distanceToOut(c7, ponr2, vmx) ≈ 100/√2 - √(95^2 - 100^2/2)
    @test distanceToOut(c8a, pparr2, vparr) ≈ 100*√5/2
    @test distanceToOut(c8a, pparr2, -vparr) ≈ 0.
    @test distanceToOut(c8a, pparr2, vz) ≈ 100.
    @test distanceToOut(c8a, pparr2, vmz) ≈ 0.
    @test distanceToOut(c8a, pparr3, vparr) ≈ 0.
    @test distanceToOut(c8a, pparr3, -vparr) ≈ 100*√5/2
    @test distanceToOut(c8a, pparr3, vz) ≈ 0.
    @test distanceToOut(c8a, pparr3, vmz) ≈ 50.
    @test distanceToOut(c8b, pparr2, vparr) ≈ 100*√5/2
    @test distanceToOut(c8b, pparr2, -vparr) ≈ 0
    @test distanceToOut(c8b, pparr2, vz) ≈ 50
    @test distanceToOut(c8b, pparr2, vmz) ≈ 0
    @test distanceToOut(c8b, pparr3, vparr) ≈ 0
    @test distanceToOut(c8b, pparr3, -vparr) ≈ 100*√5/2
    @test distanceToOut(c8b, pparr3, vz) ≈ 0
    @test distanceToOut(c8b, pparr3, vmz) ≈ 50
    @test distanceToOut(c9, Point3{T}(1e3 * kTolerance(T), 0, 50), vx2mz) ≈ 111.8033988
    @test distanceToOut(c9, Point3{T}(5, 0, 50), vx2mz) ≈ 111.8033988
    @test distanceToOut(c9, Point3{T}(10, 0, 50), vx2mz) ≈ 111.8033988
    point = Point3{T}(0.28628920024909, -0.43438111004815, -2949.0)
    dir   = Vector3{T}(6.0886686196674e-05, -9.2382200635766e-05, 0.99999999387917)
    tol   = Vector3{T}(0, 0, 0.25e3 * kTolerance(T))
    @test distanceToOut(cms, point, dir) ≈ 5898.0
    @test distanceToOut(cms, point+tol, dir) ≈ 5898.0
    #@test distanceToOut(cms, point-tol, dir) ≈ 5898.0
    point = Point3{T}(-344.13684353113, 258.98049377272, -158.20772167926)
    dir   = Vector3{T}(-0.30372024336672, -0.5581146924652, 0.77218003329776)
    @test distanceToOut(cms2, point, dir) ≈ 0.
    #@test distanceToOut(ctest10, pct10e2, d1) ≈ 111.8033988
    #@test distanceToOut(ctest10, pct10e3, d1) ≈ 111.8033988

    # distanceToIn
    @test distanceToIn(c1, pplz, vmz) ≈ Inf
    @test distanceToIn(c2, pplz, vmz) ≈ Inf
    @test distanceToIn(c3, pplz, vmz) ≈ Inf
    @test distanceToIn(c4, pplz, vmz) ≈ Inf
    @test distanceToIn(c5, pplz, vmz) ≈ Inf
    @test distanceToIn(c6, pplz, vmz) ≈ 70.
    @test distanceToIn(c7, pplz, vmz) ≈ Inf
    @test distanceToIn(c8a, pplz, vmz) ≈ 70.
    @test distanceToIn(c8b, pplz, vmz) ≈ Inf
    @test distanceToIn(c8c, pplz, vmz) ≈ Inf
    @test distanceToIn(c9, pplz, vmz) ≈ 70.
    @test distanceToIn(c9, Point3{T}(0,0,50), vmz) ≈ 0.

    @test distanceToIn(c1, pmiz, vz) ≈ Inf
    @test distanceToIn(c2, pmiz, vz) ≈ Inf
    @test distanceToIn(c3, pmiz, vz) ≈ Inf
    @test distanceToIn(c4, pmiz, vz) ≈ Inf
    @test distanceToIn(c5, pmiz, vz) ≈ Inf
    @test distanceToIn(c6, pmiz, vz) ≈ 70.
    @test distanceToIn(c7, pmiz, vz) ≈ Inf
    @test distanceToIn(c8a, pmiz, vz) ≈ 70.
    @test distanceToIn(c8b, pmiz, vz) ≈ Inf
    @test distanceToIn(c8c, pmiz, vz) ≈ Inf
    @test distanceToIn(c9, pmiz, vz) ≈ Inf
    @test distanceToIn(c1, pplx, vmx) ≈ 20
    @test distanceToIn(c1, pplz, vx) ≈ Inf
    @test distanceToIn(c4, pply, vmy) ≈ Inf
    @test distanceToIn(c1, pydx, vmy) ≈ 70
    @test distanceToIn(c3, pydx, vmy) ≈ 150 - 60tan(π/6)
    @test distanceToIn(c1, pplx, vmx) ≈ 20
    @test distanceToIn(c1, pplx, vx) ≈ Inf
    @test distanceToIn(c4, pbigx, vmx) ≈ 350
    @test distanceToIn(c4, pzero, vx) ≈ 50
    @test distanceToIn(c1, ponr2, vx) ≈ Inf
    @test distanceToIn(c1, ponr2, vmx) ≈ 0
    @test distanceToIn(c1, ponr1, vx) ≈ 0
    @test distanceToIn(c1, ponr1, vmx) ≈ 100/√2
    @test distanceToIn(c1, ponr2, vmxmy) ≈ 0
    @test distanceToIn(c8b, pparr1, vparr) ≈ 100√5/2
    @test distanceToIn(c8b, pparr2, -vparr) ≈ Inf
    @test distanceToIn(c8b, pparr3, vparr) ≈ Inf
    @test distanceToIn(c8b, pparr3, -vparr) ≈ 0
    @test distanceToIn(c8a, proot1, vz) ≈ 1000
    @test distanceToIn(c8c, proot2, vz) ≈ 1000
    #@test distanceToIn(cms2, Point3{T}(-344.13684353113, 258.98049377272, -158.20772167926), 
    #                        Vector3{T}(-0.30372022869765, -0.55811472925794, 0.77218001247454)) ≈ Inf
    @test distanceToIn(ctest10, pct10, vx) ≈ Inf
    @test distanceToIn(ctest10, pct10, vmx) ≈ 110
    @test distanceToIn(ctest10, pct10, vy) ≈ 10.579618842507896
    @test distanceToIn(ctest10, pct10, vmy) ≈ 71.50521555565261
    @test distanceToIn(ctest10, pct10, vz) ≈ Inf
    @test distanceToIn(ctest10, pct10phi1, vx) ≈ Inf
    @test isapprox(distanceToIn(ctest10, pct10phi1, vmx), 0., atol=100*eps(T)) 
    @test isapprox(distanceToIn(ctest10, pct10phi1, vy), 0., atol=100*eps(T))
    @test distanceToIn(ctest10, pct10phi1, vmy) ≈ 80.83778132003165
    @test distanceToIn(ctest10, pct10phi1, vz) ≈ Inf
    @test distanceToIn(ctest10, pct10phi1, vmz) ≈ Inf
    @test distanceToIn(ctest10, pct10phi2, vx) ≈ Inf
    @test isapprox(distanceToIn(ctest10, pct10phi2, vmx), 0., atol=100*eps(T))
    @test distanceToIn(ctest10, pct10phi2, vy) ≈ 77.78352
    @test isapprox(distanceToIn(ctest10, pct10phi2, vmy), 0., atol=100*eps(T))
    @test distanceToIn(ctest10, pct10phi2, vz) ≈ Inf
    @test distanceToIn(ctest10, pct10phi2, vmz) ≈ Inf
    @test distanceToIn(ctest10, pct10mx, vx) ≈ Inf
    @test distanceToIn(ctest10, pct10mx, vmx) ≈ 0.
    @test distanceToIn(ctest10, pct10mx, vy) ≈ Inf
    @test distanceToIn(ctest10, pct10mx, vmy) ≈ Inf
    @test distanceToIn(ctest10, pct10mx, vz) ≈ Inf
    @test distanceToIn(ctest10, pct10mx, vmz) ≈ 0.
    @test distanceToIn(ctest10, pct10e1, d1) ≈ 1171.1712376073724
    @test distanceToIn(ctest10, pct10e4, d1) ≈ Inf
    @test distanceToIn(ctest10, pt10s2, vt10d) ≈ 142.80171070987436
    @test distanceToIn(ctest10, Point3{T}(90cos(10deg), 90sin(10deg), 0), vmy) ≈ Inf
    @test distanceToIn(ctest10, Point3{T}(90cos(50deg),-90sin(50deg), 0), vx) ≈ Inf
    alex1P = Point3{T}(49.840299921054168, -59.39735648688918, -20.893051766050633)
    alex1V = Vector3{T}(0.6068108874999103, 0.35615926907657169, 0.71058505603651234)
    @test distanceToIn(ctest10, alex1P, alex1V) ≈ Inf
    @test distanceToOut(ctest10, alex1P, alex1V) ≈ 0
    alex2P = Point3{T}(127.0075852717127, -514.1050841937065, 69.47104834263656)
    alex2V = Vector3{T}(0.1277616879490939, 0.4093610465777845, 0.9033828007202369)
    @test inside(ctest10, alex2P) == kOutside
    @test distanceToIn(ctest10, alex2P, alex2V) ≈ Inf
end
