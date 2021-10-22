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

    #@test distanceToOut(trd1, Point3{T}(0,0,0), vx) ≈ 20


end

#=
 // DistanceToOut(P,V)

  Dist = trd1.DistanceToOut(pzero, vx);
  assert(ApproxEqual(Dist, 20));
  Dist = trd1.DistanceToOut(pzero, vmx);
  assert(ApproxEqual(Dist, 20));
  Dist = trd1.DistanceToOut(pzero, vy);
  assert(ApproxEqual(Dist, 30));
  Dist = trd1.DistanceToOut(pzero, vmy);
  assert(ApproxEqual(Dist, 30));
  Dist = trd1.DistanceToOut(pzero, vz);
  assert(ApproxEqual(Dist, 40));
  Dist = trd1.DistanceToOut(pzero, vmz);
  assert(ApproxEqual(Dist, 40));
  Dist = trd1.DistanceToOut(pzero, vxy);
  assert(ApproxEqual(Dist, std::sqrt(800.)));

  Dist = trd1.DistanceToOut(ponxside, vx);
  assert(ApproxEqual(Dist, 0));
  Dist = trd1.DistanceToOut(ponmxside, vmx);
  assert(ApproxEqual(Dist, 0));
  Dist = trd1.DistanceToOut(ponyside, vy);
  assert(ApproxEqual(Dist, 0));
  Dist = trd1.DistanceToOut(ponmyside, vmy);
  assert(ApproxEqual(Dist, 0));
  Dist = trd1.DistanceToOut(ponzside, vz);
  assert(ApproxEqual(Dist, 0));
  Dist = trd1.DistanceToOut(ponmzside, vmz);
  assert(ApproxEqual(Dist, 0));

  Dist  = trd2.DistanceToOut(pzero, vx);
  valid = trd2.Normal(pzero + Dist * vx, normal);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(normal, Vec_t(cosa, 0, -sina)));

  Dist  = trd2.DistanceToOut(pzero, vmx);
  valid = trd2.Normal(pzero + Dist * vmx, normal);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(normal, Vec_t(-cosa, 0, -sina)));

  Dist  = trd2.DistanceToOut(pzero, vy);
  valid = trd2.Normal(pzero + Dist * vy, normal);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(normal, Vec_t(0, cosa, -sina)));

  Dist  = trd2.DistanceToOut(pzero, vmy);
  valid = trd2.Normal(pzero + Dist * vmy, normal);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(normal, Vec_t(0, -cosa, -sina)));

  Dist  = trd2.DistanceToOut(pzero, vz);
  valid = trd2.Normal(pzero + Dist * vz, normal);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(normal, vz));

  Dist  = trd2.DistanceToOut(pzero, vmz);
  valid = trd2.Normal(pzero + Dist * vmz, normal);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(normal, vmz));

  Dist  = trd2.DistanceToOut(pzero, vxy);
  valid = trd2.Normal(pzero + Dist * vxy, normal);
  assert(ApproxEqual(Dist, std::sqrt(800.)));

  Dist  = trd2.DistanceToOut(ponxside, vx);
  valid = trd2.Normal(ponxside + Dist * vx, normal);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(normal, Vec_t(cosa, 0, -sina)));

  Dist  = trd2.DistanceToOut(ponmxside, vmx);
  valid = trd2.Normal(ponmxside + Dist * vmx, normal);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(normal, Vec_t(-cosa, 0, -sina)));

  Dist  = trd2.DistanceToOut(ponyside, vy);
  valid = trd2.Normal(ponyside + Dist * vy, normal);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(normal, Vec_t(0, cosa, -sina)));

  Dist  = trd2.DistanceToOut(ponmyside, vmy);
  valid = trd2.Normal(ponmyside + Dist * vmy, normal);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(normal, Vec_t(0, -cosa, -sina)));

  Dist  = trd2.DistanceToOut(ponzside, vz);
  valid = trd2.Normal(ponzside + Dist * vz, normal);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(normal, vz));

  Dist  = trd2.DistanceToOut(ponmzside, vmz);
  valid = trd2.Normal(ponmzside + Dist * vmz, normal);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(normal, vmz));
=#