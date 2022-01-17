@testset "Tube{$T}" for T in (Float64, Float32)
    pzero = Point3{T}(0,0,0)
    pbigx = Point3{T}(100,0,0)
    pbigy = Point3{T}(0,100,0)
    pbigz = Point3{T}(0,0,100)
    pbigmx = Point3{T}(-100,0,0)
    pbigmy = Point3{T}(0,-100,0)
    pbigmz = Point3{T}(0,0,-100)
    ponxside = Point3{T}(50,0,0)
    ponyside = Point3{T}(0,50,0)
    ponzside = Point3{T}(0,0,50)
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
    
    t1  = Tube{T}(0, 50, 50, 0, 2π)
    t1a = Tube{T}(0, 50, 50, 0, π/2)
    t2  = Tube{T}(45, 50, 50, 0, 2π)
    t2a = Tube{T}(5, 50, 50, 0, 2π)
    t2b = Tube{T}(15, 50, 50, 0, 2π)
    t2c = Tube{T}(25, 50, 50, 0, 2π)
    t2d = Tube{T}(35, 50, 50, 0, 2π)
    t3  = Tube{T}(0, 50, 50, π/2, π/2)
    t4  = Tube{T}(45, 50, 50, π/2, π/2)
    t5  = Tube{T}(50, 100, 50, 0.0, 3π/2)
    t6  = Tube{T}(0, 50, 50, π/2, 3π/2);
    tube6 = Tube{T}(750, 760, 350, 0.31415926535897931, 5.6548667764616276)
    tube7 = Tube{T}(2200, 3200, 2500, -0.68977164349384879, 3.831364227270472)
    tube8 = Tube{T}(2550, 2580, 2000, 0, 2π)
    tube9 = Tube{T}(1150, 1180, 2000, 0, 2π)
    tube10 = Tube{T}(400, 405, 400, 0, 2π)
    clad = Tube{T}(90., 110., 105, 0., π)
    core = Tube{T}(95., 105., 100, 0., π)

    @test capacity(t1) ≈ 50 * 2 * π * 50 * 50
    @test surface(t2) ≈ T(2) * π * (45 + 50) * (50 - 45 + 2 * 50)

    # Check Inside
    @test inside(t1, pzero) == kInside
    @test inside(t1, pbigz) == kOutside
    @test inside(t1, ponxside) == kSurface
    @test inside(t1, ponyside) == kSurface
    @test inside(t1, ponzside) == kSurface
    @test inside(t1a, pzero) == kSurface
    @test inside(t1a, pbigz) == kOutside
    @test inside(t1a, ponxside) == kSurface
    @test inside(t1a, ponyside) == kSurface
    @test inside(t1a, ponzside) == kSurface

    @test inside(t2, pzero) == kOutside
    @test inside(t2, pbigz) == kOutside
    @test inside(t2, ponxside) == kSurface
    @test inside(t2, ponyside) == kSurface
    @test inside(t2, ponzside) == kOutside
    @test inside(t2a, pzero) == kOutside
    @test inside(t2a, pbigz) == kOutside
    @test inside(t2a, ponxside) == kSurface
    @test inside(t2a, ponyside) == kSurface
    @test inside(t2a, ponzside) == kOutside

    # Check Normal

    @test normal(t1, ponxside) ≈ vx
    @test normal(t4, Point3{T}(0,50,0)) ≈ Vector3{T}(1/√2, 1/√2, 0)
    @test normal(t4, Point3{T}(0,45,0)) ≈ Vector3{T}(1/√2, -1/√2, 0)
    @test normal(t4, Point3{T}(0,45,50)) ≈ Vector3{T}(1/√3, -1/√3, 1/√3)
    @test normal(t4, Point3{T}(0,45,-50)) ≈ Vector3{T}(1/√3, -1/√3, -1/√3)
    @test normal(t4, Point3{T}(-50,0,-50)) ≈ Vector3{T}(-1/√3, -1/√3, -1/√3)
    @test normal(t4, Point3{T}(-50,0,0)) ≈ Vector3{T}(-1/√2, -1/√2, 0)
    @test normal(t6, Point3{T}(0,0,0)) ≈ Vector3{T}(1/√2, 1/√2, 0)

    # SafetyToOut(P)
    @test safetyToOut(t1, pzero) ≈ 50

    # DistanceToOut(P,V)
    dist = distanceToOut(t1, pzero, vx)
    @test dist ≈ 50 && normal(t1, pzero + dist * vx) ≈ vx

    dist = distanceToOut(t1, pzero, vmx)
    @test dist ≈ 50 && normal(t1, pzero + dist * vmx) ≈ vmx

    dist = distanceToOut(t1, pzero, vy)
    @test dist ≈ 50 && normal(t1, pzero + dist * vy) ≈ vy

    dist = distanceToOut(t1, pzero, vmy)
    @test dist ≈ 50 && normal(t1, pzero + dist * vmy) ≈ vmy

    dist = distanceToOut(t1, pzero, vz)
    @test dist ≈ 50 && normal(t1, pzero + dist * vz) ≈ vz

    dist = distanceToOut(t1, pzero, vmz)
    @test dist ≈ 50 && normal(t1, pzero + dist * vmz) ≈ vmz

    dist = distanceToOut(t1, pzero, vxy)
    @test dist ≈ 50 && normal(t1, pzero + dist * vxy) ≈ vxy

    @test isapprox(distanceToOut(t3, Point3{T}(0, 10, 0), vx), 0., atol=10*eps(T))
    @test distanceToOut(t3, Point3{T}(-0.5, 10, 0), vx) ≈ 0.5
    @test distanceToOut(t3, Point3{T}(-0.5, 9, 0), vx) ≈ 0.5
    @test distanceToOut(t3, Point3{T}(-5, 9.5, 0), vx) ≈ 5.
    @test distanceToOut(t3, Point3{T}(-5, 9.5, 0), vmx) ≈ 44.089204515860
    @test distanceToOut(t3, Point3{T}(-5, 9, 0), vxmy) ≈ 7.0710678

    # SafetyToIn(P)
    @test safetyToIn(t1, pbigx) ≈ 50
    @test safetyToIn(t1, pbigmx) ≈ 50
    @test safetyToIn(t1, pbigy) ≈ 50
    @test safetyToIn(t1, pbigmy) ≈ 50
    @test safetyToIn(t1, pbigz) ≈ 50
    @test safetyToIn(t1, pbigmz) ≈ 50

    # DistanceToIn(P,V)
    @test distanceToIn(t1, pbigx, vmx) ≈ 50
    @test distanceToIn(t1, pbigmx, vx) ≈ 50
    @test distanceToIn(t1, pbigy, vmy) ≈ 50
    @test distanceToIn(t1, pbigmy, vy) ≈ 50
    @test distanceToIn(t1, pbigz, vmz) ≈ 50
    @test distanceToIn(t1, pbigmz, vz) ≈ 50
    @test distanceToIn(t1, pbigx, vxy) == Inf

    @test distanceToIn(t1a, pbigz, vmz) ≈ 50
    @test distanceToOut(t2, Point3{T}(45.5, 0, 0), vx) ≈ 4.5
    @test distanceToOut(t2, Point3{T}(45.5, 0, 0), vmx) ≈ 0.5
    @test distanceToOut(t2, Point3{T}(49.5, 0, 0), vmx) ≈ 4.5
    @test distanceToOut(t2, Point3{T}(49.5, 0, 0), vx) ≈ 0.5

end
