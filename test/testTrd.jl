@testset "Trd{$T}" for T in (Float64, Float32)

    # Define the nededed trd's    
    trd1 = Trd{T}(20, 20, 30, 30, 40)
    trd2 = Trd{T}(10, 30, 20, 40, 40)
    trd3 = Trd{T}( 0.14999999999999999, 0.14999999999999999, 24.707000000000001, 24.707000000000001,
    22.699999999999999)

    # Volume and Surfaces
    @test capacity(trd1) ≈ 8. * 20. * 30. * 40.
    @test surface(trd1) ≈ 20800.

    # Inside
    @test inside(trd1, Point3{T}(0, 0, 0)) == kInside
    @test inside(trd1, Point3{T}(0, 0, 100)) == kOutside
    @test inside(trd1, Point3{T}(20, 0, 0)) == kSurface
    @test inside(trd1, Point3{T}(0, 30, 0)) == kSurface
    @test inside(trd1, Point3{T}(0, 0, 40)) == kSurface

    @test inside(trd1, Point3{T}(20, 30, 40)) == kSurface
    @test inside(trd1, Point3{T}(-20, 30, 40)) == kSurface
    @test inside(trd1, Point3{T}(20, -30, 40)) == kSurface
    @test inside(trd1, Point3{T}(20, 30, -40)) == kSurface
    @test inside(trd1, Point3{T}(20, 30,  0)) == kSurface
    @test inside(trd1, Point3{T}( 0, 30, 40)) == kSurface
    @test inside(trd1, Point3{T}(20,  0, 40)) == kSurface
    @test inside(trd1, Point3{T}(-20, -30, -40)) == kSurface
    
    @test inside(trd2, Point3{T}(0, 0, 0)) == kInside
    @test inside(trd2, Point3{T}(0, 0, 100)) == kOutside
    @test inside(trd2, Point3{T}(20, 0, 0)) == kSurface
    @test inside(trd2, Point3{T}(0, 30, 0)) == kSurface
    @test inside(trd2, Point3{T}(0, 0, 40)) == kSurface

    # Surface Normal
    @test normal(trd1, Point3{T}(20, 0, 0)) ≈ Vector3{T}(1,0,0)
    @test normal(trd1, Point3{T}(-20, 0, 0)) ≈ Vector3{T}(-1,0,0)
    @test normal(trd1, Point3{T}(0, 30, 0)) ≈ Vector3{T}(0,1,0)
    @test normal(trd1, Point3{T}(0, -30, 0)) ≈ Vector3{T}(0,-1,0)
    @test normal(trd1, Point3{T}(0, 0, 40)) ≈ Vector3{T}(0,0,1)
    @test normal(trd1, Point3{T}(0, 0, -40)) ≈ Vector3{T}(0,0,-1)
    @test normal(trd1, Point3{T}(0, 25,  40)) ≈ Vector3{T}(0,0,1)
    @test normal(trd1, Point3{T}(0, 25, -40)) ≈ Vector3{T}(0,0,-1)

    # Normals on Edges
    @test normal(trd1, Point3{T}( 20,  30, 0)) ≈ Vector3{T}( 1/√2,  1/√2, 0)
    @test normal(trd1, Point3{T}(-20, -30, 0)) ≈ Vector3{T}(-1/√2, -1/√2, 0)
    @test normal(trd1, Point3{T}( 20, -30, 0)) ≈ Vector3{T}( 1/√2, -1/√2, 0)
    @test normal(trd1, Point3{T}(-20,  30, 0)) ≈ Vector3{T}(-1/√2,  1/√2, 0)
    @test normal(trd1, Point3{T}( 20, 0,  40)) ≈ Vector3{T}( 1/√2, 0,  1/√2)
    @test normal(trd1, Point3{T}(-20, 0, -40)) ≈ Vector3{T}(-1/√2, 0, -1/√2)
    @test normal(trd1, Point3{T}( 20, 0, -40)) ≈ Vector3{T}( 1/√2, 0, -1/√2)
    @test normal(trd1, Point3{T}(-20, 0,  40)) ≈ Vector3{T}(-1/√2, 0,  1/√2)
    @test normal(trd1, Point3{T}(0,  30,  40)) ≈ Vector3{T}(0,  1/√2,  1/√2)
    @test normal(trd1, Point3{T}(0, -30, -40)) ≈ Vector3{T}(0, -1/√2, -1/√2)
    @test normal(trd1, Point3{T}(0, -30,  40)) ≈ Vector3{T}(0, -1/√2,  1/√2)
    @test normal(trd1, Point3{T}(0,  30, -40)) ≈ Vector3{T}(0,  1/√2, -1/√2)

    # Normals on Corners
    @test normal(trd1, Point3{T}( 20,  30,  40)) ≈ Vector3{T}( 1/√3,  1/√3,  1/√3)
    @test normal(trd1, Point3{T}(-20,  30,  40)) ≈ Vector3{T}(-1/√3,  1/√3,  1/√3)
    @test normal(trd1, Point3{T}( 20, -30,  40)) ≈ Vector3{T}( 1/√3, -1/√3,  1/√3)
    @test normal(trd1, Point3{T}(-20, -30,  40)) ≈ Vector3{T}(-1/√3, -1/√3,  1/√3)
    @test normal(trd1, Point3{T}( 20,  30, -40)) ≈ Vector3{T}( 1/√3,  1/√3, -1/√3)
    @test normal(trd1, Point3{T}(-20,  30, -40)) ≈ Vector3{T}(-1/√3,  1/√3, -1/√3)
    @test normal(trd1, Point3{T}( 20, -30, -40)) ≈ Vector3{T}( 1/√3, -1/√3, -1/√3)
    @test normal(trd1, Point3{T}(-20, -30, -40)) ≈ Vector3{T}(-1/√3, -1/√3, -1/√3)

    # Normals for trd2
    @test normal(trd2, Point3{T}( 20, 0, 0)) ≈ Vector3{T}( 4/√17, 0, -1/√17)
    @test normal(trd2, Point3{T}(-20, 0, 0)) ≈ Vector3{T}(-4/√17, 0, -1/√17)
    @test normal(trd2, Point3{T}( 0, 30, 0)) ≈ Vector3{T}(0,  4/√17, -1/√17)
    @test normal(trd2, Point3{T}( 0,-30, 0)) ≈ Vector3{T}(0, -4/√17, -1/√17)
    @test normal(trd2, Point3{T}( 0, 0, 40)) ≈ Vector3{T}(0,  0,  1)
    @test normal(trd2, Point3{T}( 0, 0,-40)) ≈ Vector3{T}(0,  0, -1)
    @test normal(trd2, Point3{T}( 0, 25, 40)) ≈ Vector3{T}(0,  0,  1)
    @test normal(trd2, Point3{T}( 0, 25,-40)) ≈ Vector3{T}(0,  0, -1)


    # SafetyToOut(P)
    @test safetyToOut(trd1, Point3{T}(0,0,0)) ≈ 20
    @test safetyToOut(trd1, Point3{T}(1,0,0)) ≈ 19
    @test safetyToOut(trd1, Point3{T}(0,1,0)) ≈ 20
    @test safetyToOut(trd1, Point3{T}(0,0,1)) ≈ 20

    @test safetyToOut(trd2, Point3{T}(0,0,0)) ≈ 20*4/√17
    @test safetyToOut(trd2, Point3{T}(1,0,0)) ≈ 19*4/√17
    @test safetyToOut(trd2, Point3{T}(0,1,0)) ≈ 20*4/√17
    @test safetyToOut(trd2, Point3{T}(0,0,1)) ≈ 20*4/√17 + 1/√17

    # SafetyToIn(P)
    @test safetyToIn(trd1, Point3{T}(100,0,0)) ≈ 80
    @test safetyToIn(trd1, Point3{T}(-100,0,0)) ≈ 80
    @test safetyToIn(trd1, Point3{T}(0, 100,0)) ≈ 70
    @test safetyToIn(trd1, Point3{T}(0,-100,0)) ≈ 70
    @test safetyToIn(trd1, Point3{T}(0,0, 100)) ≈ 60
    @test safetyToIn(trd1, Point3{T}(0,0,-100)) ≈ 60

    @test safetyToIn(trd2, Point3{T}( 100,0,0)) ≈ 80 * 4/√17
    @test safetyToIn(trd2, Point3{T}(-100,0,0)) ≈ 80 * 4/√17
    @test safetyToIn(trd2, Point3{T}(0, 100,0)) ≈ 70 * 4/√17
    @test safetyToIn(trd2, Point3{T}(0,-100,0)) ≈ 70 * 4/√17
    @test safetyToIn(trd2, Point3{T}(0,0, 100)) ≈ 60
    @test safetyToIn(trd2, Point3{T}(0,0,-100)) ≈ 60
  
    # DistanceToOut(P,V)

    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}( 1, 0, 0)) ≈ 20
    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}(-1, 0, 0)) ≈ 20
    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}(0,  1, 0)) ≈ 30
    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}(0, -1, 0)) ≈ 30
    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}(0, 0,  1)) ≈ 40
    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}(0, 0, -1)) ≈ 40
    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}(0,  1, 0)) ≈ 30
    @test distanceToOut(trd1, Point3{T}(0,0,0), Vector3{T}(1/√2, 1/√2, 0.)) ≈ √T(800)

    @test distanceToOut(trd1, Point3{T}( 20,0,0), Vector3{T}( 1, 0, 0)) ≈ 0
    @test distanceToOut(trd1, Point3{T}(-20,0,0), Vector3{T}(-1, 0, 0)) ≈ 0
    @test distanceToOut(trd1, Point3{T}(0, 30,0), Vector3{T}( 0, 1, 0)) ≈ 0
    @test distanceToOut(trd1, Point3{T}(0,-30,0), Vector3{T}( 0,-1, 0)) ≈ 0
    @test distanceToOut(trd1, Point3{T}(0,0, 40), Vector3{T}( 0, 0, 1)) ≈ 0
    @test distanceToOut(trd1, Point3{T}(0,0,-40), Vector3{T}( 0, 0,-1)) ≈ 0

    dist = distanceToOut(trd2, Point3{T}(0,0,0), Vector3{T}( 1, 0, 0))
    @test dist ≈ 20 && normal(trd2, Point3{T}(0,0,0) + Vector3{T}( 1, 0, 0) * dist) ≈ Vector3{T}( 4/√17, 0, -1/√17)
    dist = distanceToOut(trd2, Point3{T}(0,0,0), Vector3{T}(-1, 0, 0))
    @test dist ≈ 20 && normal(trd2, Point3{T}(0,0,0) + Vector3{T}(-1, 0, 0) * dist) ≈ Vector3{T}(-4/√17, 0, -1/√17)
 
    dist = distanceToOut(trd2, Point3{T}(0,0,0), Vector3{T}( 0, 1, 0))
    @test dist ≈ 30 && normal(trd2, Point3{T}(0,0,0) + Vector3{T}( 0, 1, 0) * dist) ≈ Vector3{T}( 0, 4/√17,-1/√17)
    dist = distanceToOut(trd2, Point3{T}(0,0,0), Vector3{T}( 0,-1, 0))
    @test dist ≈ 30 && normal(trd2, Point3{T}(0,0,0) + Vector3{T}( 0,-1, 0) * dist) ≈ Vector3{T}( 0,-4/√17,-1/√17)

    dist = distanceToOut(trd2, Point3{T}(0,0,0), Vector3{T}( 0, 0, 1))
    @test dist ≈ 40 && normal(trd2, Point3{T}(0,0,0) + Vector3{T}( 0, 0, 1) * dist) ≈ Vector3{T}( 0, 0, 1)
    dist = distanceToOut(trd2, Point3{T}(0,0,0), Vector3{T}( 0, 0,-1))
    @test dist ≈ 40 && normal(trd2, Point3{T}(0,0,0) + Vector3{T}( 0, 0,-1) * dist) ≈ Vector3{T}( 0, 0,-1)

    @test dist = distanceToOut(trd2, Point3{T}(0,0,0), Vector3{T}(1/√2, 1/√2, 0.)) ≈  √T(800)

    dist = distanceToOut(trd2, Point3{T}( 20,0,0), Vector3{T}( 1, 0, 0))
    @test dist ≈ 0 && normal(trd2, Point3{T}( 20,0,0) + Vector3{T}( 1, 0, 0) * dist) ≈ Vector3{T}( 4/√17, 0, -1/√17)
    dist = distanceToOut(trd2, Point3{T}(-20,0,0), Vector3{T}(-1, 0, 0))
    @test dist ≈ 0 && normal(trd2, Point3{T}(-20,0,0) + Vector3{T}(-1, 0, 0) * dist) ≈ Vector3{T}(-4/√17, 0, -1/√17)

    dist = distanceToOut(trd2, Point3{T}( 0, 30, 0), Vector3{T}( 0, 1, 0))
    @test dist ≈ 0 && normal(trd2, Point3{T}( 0, 30, 0) + Vector3{T}( 0, 1, 0) * dist) ≈ Vector3{T}( 0, 4/√17, -1/√17)
    dist = distanceToOut(trd2, Point3{T}(0,-30,0), Vector3{T}(0, -1, 0))
    @test dist ≈ 0 && normal(trd2, Point3{T}(0,-30,0) + Vector3{T}( 0, -1, 0) * dist) ≈ Vector3{T}( 0, -4/√17, -1/√17)

    dist = distanceToOut(trd2, Point3{T}( 0, 0, 40), Vector3{T}( 0, 0, 1))
    @test dist ≈ 0 && normal(trd2, Point3{T}( 0, 0, 40) + Vector3{T}( 0, 0, 1) * dist) ≈ Vector3{T}( 0, 0, 1)
    dist = distanceToOut(trd2, Point3{T}(0, 0,-40), Vector3{T}(0, 0,-1))
    @test dist ≈ 0 && normal(trd2, Point3{T}( 0, 0,-40) + Vector3{T}( 0, 0, -1) * dist) ≈ Vector3{T}( 0, 0, -1)

    # DistanceToIn(P,V)

    @test distanceToIn(trd1, Point3{T}( 100, 0, 0), Vector3{T}(-1, 0, 0)) ≈ 80
    @test distanceToIn(trd1, Point3{T}(-100, 0, 0), Vector3{T}( 1, 0, 0)) ≈ 80
    @test distanceToIn(trd1, Point3{T}( 0, 100, 0), Vector3{T}( 0,-1, 0)) ≈ 70
    @test distanceToIn(trd1, Point3{T}( 0,-100, 0), Vector3{T}( 0, 1, 0)) ≈ 70
    @test distanceToIn(trd1, Point3{T}( 0, 0, 100), Vector3{T}( 0, 0,-1)) ≈ 60
    @test distanceToIn(trd1, Point3{T}( 0, 0,-100), Vector3{T}( 0, 0, 1)) ≈ 60
    @test distanceToIn(trd1, Point3{T}( 100, 0, 0), Vector3{T}( 1/√2, 1/√2, 0)) ≈ Inf
    @test distanceToIn(trd1, Point3{T}(-100, 0, 0), Vector3{T}( 1/√2, 1/√2, 0)) ≈ Inf

    @test distanceToIn(trd2, Point3{T}( 100, 0, 0), Vector3{T}(-1, 0, 0)) ≈ 80
    @test distanceToIn(trd2, Point3{T}(-100, 0, 0), Vector3{T}( 1, 0, 0)) ≈ 80
    @test distanceToIn(trd2, Point3{T}( 0, 100, 0), Vector3{T}( 0,-1, 0)) ≈ 70
    @test distanceToIn(trd2, Point3{T}( 0,-100, 0), Vector3{T}( 0, 1, 0)) ≈ 70
    @test distanceToIn(trd2, Point3{T}( 0, 0, 100), Vector3{T}( 0, 0,-1)) ≈ 60
    @test distanceToIn(trd2, Point3{T}( 0, 0,-100), Vector3{T}( 0, 0, 1)) ≈ 60
    @test distanceToIn(trd2, Point3{T}( 100, 0, 0), Vector3{T}( 1/√2, 1/√2, 0)) ≈ Inf
    @test distanceToIn(trd2, Point3{T}(-100, 0, 0), Vector3{T}( 1/√2, 1/√2, 0)) ≈ Inf


    tdir = normalize(Vector3{T}(-0.76165597579890043, 0.64364445891356026, -0.074515708658524193))
    #@test distanceToIn(trd3, Point3{T}( 0.15000000000000185, -22.048743592955137, 2.4268539333219472), tdir) ≈ 0
      
end
