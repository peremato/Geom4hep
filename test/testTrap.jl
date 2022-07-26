@testset "Trap{$T}" for T in (Float64, Float32)
    trap1 = Trap{T}(40, 0, 0, 30, 20, 20, 0, 30, 20, 20, 0) #  box:20,30,40
    trap2 = Trap{T}(40, 0, 0, 20, 10, 10, 0, 40, 30, 30, 0)
    trap3 = Trap{T}(50, 0, 0, 50, 50, 50, π/4, 50, 50, 50, π/4)
    trap4 = Trap{T}(50, 0, 0, 50, 50, 50, -π/4, 50, 50, 50, -π/4)
    trap5 = Trap{T}(20, 20, 30, 30, 40)  # trd like construction

    # Check cubic volume
    @test capacity(trap1) ≈ 8 * 20 * 30 * 40
    @test capacity(trap2) ≈ 2 * 40. * ((20. + 40.) * (10. + 30.) + (30. - 10.) * (40. - 20.) / 3.)
    @test capacity(trap3) ≈ 8 * 50 * 50 * 50
    @test capacity(trap4) ≈ 8 * 50 * 50 * 50

    # Check surface
    @test surface(trap1) ≈ 2 * (40 * 60 + 80 * 60 + 80 * 40)
    @test surface(trap2) ≈ 4 * (10 * 20 + 30 * 40) + 2 * ((20 + 40) * sqrt(4 * 40 * 40 + (30 - 10) * (30 - 10)) +
                                                          (30 + 10) * sqrt(4 * 40 * 40 + (40 - 20) * (40 - 20)))
    @test surface(trap5) ≈ 20800

    pzero = Point3{T}(0,0,0)
    pbig  = Point3{T}(100, 100, 100)
    pbigx = Point3{T}(100,0,0)
    pbigy = Point3{T}(0,100,0)
    pbigz = Point3{T}(0,0,100)
    pbigmx = Point3{T}(-100,0,0)
    pbigmy = Point3{T}(0,-100,0)
    pbigmz = Point3{T}(0,0,-100)
    ponxside = Point3{T}(20,0,0)
    ponyside = Point3{T}(0,30,0)
    ponzside = Point3{T}(0,0,40)
    ponmxside = Point3{T}(-20,0,0)
    ponmyside = Point3{T}(0,-30,0)
    ponmzside = Point3{T}(0,0,-40)
    vx    = Vector3{T}( 1, 0, 0)
    vy    = Vector3{T}( 0, 1, 0)
    vz    = Vector3{T}( 0, 0, 1)
    vmx   = Vector3{T}(-1, 0, 0)    
    vmy   = Vector3{T}( 0,-1, 0)
    vmz   = Vector3{T}( 0, 0,-1)
    vxy   = Vector3{T}( 1/√2, 1/√2, 0)
    vmxy  = Vector3{T}(-1/√2, 1/√2, 0)
    vxmy  = Vector3{T}( 1/√2,-1/√2, 0)
    vmxmy = Vector3{T}(-1/√2,-1/√2, 0)
    vxz   = Vector3{T}( 1/√2, 0, 1/√2)
    vxmz  = Vector3{T}( 1/√2, 0,-1/√2)
    vyz   = Vector3{T}( 0, 1/√2, 1/√2)
    vymz  = Vector3{T}( 0, 1/√2, -1/√2)
    # SafetyToOut(P)
    @test safetyToOut(trap1, pzero) ≈ 20
    @test safetyToOut(trap1, Point3{T}(1,0,0)) ≈ 19
    @test safetyToOut(trap1, Point3{T}(0,1,0)) ≈ 20
    @test safetyToOut(trap1, Point3{T}(0,0,1)) ≈ 20

    @test safetyToOut(trap2, pzero) ≈ 20*4/√17
    @test safetyToOut(trap2, Point3{T}(1,0,0)) ≈ 19*4/√17
    @test safetyToOut(trap2, Point3{T}(0,1,0)) ≈ 20*4/√17
    @test safetyToOut(trap2, Point3{T}(0,0,1)) ≈ 20*4/√17 + 1/√17

    # SafetyToIn(P)
    @test safetyToIn(trap1, pbig) ≈ 80
    @test safetyToIn(trap1, pbigx) ≈ 80
    @test safetyToIn(trap1, pbigmx) ≈ 80
    @test safetyToIn(trap1, pbigy) ≈ 70
    @test safetyToIn(trap1, pbigmy) ≈ 70
    @test safetyToIn(trap1, pbigz) ≈ 60
    @test safetyToIn(trap1, pbigmz) ≈ 60

    @test safetyToIn(trap2, pbigx) ≈ 80*4/√17
    @test safetyToIn(trap2, pbigmx) ≈ 80*4/√17
    @test safetyToIn(trap2, pbigy) ≈ 70*4/√17
    @test safetyToIn(trap2, pbigmy) ≈ 70*4/√17
    @test safetyToIn(trap2, pbigz) ≈ 60
    @test safetyToIn(trap2, pbigmz) ≈ 60


    # DistanceToOut(P,V)
    @test distanceToOut(trap1, pzero, vx) ≈ 20
    @test distanceToOut(trap1, pzero, vmx) ≈ 20
    @test distanceToOut(trap1, pzero, vy) ≈ 30
    @test distanceToOut(trap1, pzero, vmy) ≈ 30
    @test distanceToOut(trap1, pzero, vz) ≈ 40
    @test distanceToOut(trap1, pzero, vmz) ≈ 40
    @test distanceToOut(trap1, pzero, vxy) ≈ √800
    @test distanceToOut(trap1, ponxside, vx) ≈ 0
    @test distanceToOut(trap1, ponmxside, vmx) ≈ 0
    @test distanceToOut(trap1, ponyside, vy) ≈ 0
    @test distanceToOut(trap1, ponmyside, vmy) ≈ 0
    @test distanceToOut(trap1, ponzside, vz) ≈ 0
    @test distanceToOut(trap1, ponmzside, vmz) ≈ 0
    @test distanceToOut(trap1, ponxside, vmx) ≈ 40
    @test distanceToOut(trap1, ponmxside, vx) ≈ 40
    @test distanceToOut(trap1, ponyside, vmy) ≈ 60
    @test distanceToOut(trap1, ponmyside, vy) ≈ 60
    @test distanceToOut(trap1, ponzside, vmz) ≈ 80
    @test distanceToOut(trap1, ponmzside, vz) ≈ 80

    @test distanceToOut(trap2, pzero, vx) ≈ 20
    @test distanceToOut(trap2, pzero, vmx) ≈ 20
    @test distanceToOut(trap2, pzero, vy) ≈ 30
    @test distanceToOut(trap2, pzero, vmy) ≈ 30
    @test distanceToOut(trap2, pzero, vz) ≈ 40
    @test distanceToOut(trap2, pzero, vmz) ≈ 40
    @test distanceToOut(trap2, pzero, vxy) ≈ √800
    @test distanceToOut(trap2, ponxside, vx) ≈ 0
    @test distanceToOut(trap2, ponmxside, vmx) ≈ 0
    @test isapprox(distanceToOut(trap2, ponyside, vy), 0, atol=1e-6)
    @test distanceToOut(trap2, ponmyside, vmy) ≈ 0
    @test distanceToOut(trap2, ponzside, vz) ≈ 0
    @test distanceToOut(trap2, ponmzside, vmz) ≈ 0

    # DistanceToIn(P,V)
    @test distanceToIn(trap1, pbigx, vmx) ≈ 80
    @test distanceToIn(trap1, pbigmx, vx) ≈ 80
    @test distanceToIn(trap1, pbigy, vmy) ≈ 70
    @test distanceToIn(trap1, pbigmy, vy) ≈ 70
    @test distanceToIn(trap1, pbigz, vmz) ≈ 60
    @test distanceToIn(trap1, pbigmz, vz) ≈ 60
    @test distanceToIn(trap1, pbigx, vxy) ≈ Inf
    @test distanceToIn(trap1, pbigmx, vxy) ≈ Inf

    @test distanceToIn(trap2, pbigx, vmx) ≈ 80
    @test distanceToIn(trap2, pbigmx, vx) ≈ 80
    @test distanceToIn(trap2, pbigy, vmy) ≈ 70
    @test distanceToIn(trap2, pbigmy, vy) ≈ 70
    @test distanceToIn(trap2, pbigz, vmz) ≈ 60
    @test distanceToIn(trap2, pbigmz, vz) ≈ 60
    @test distanceToIn(trap2, pbigx, vxy) ≈ Inf
    @test distanceToIn(trap2, pbigmx, vxy) ≈ Inf

    @test distanceToIn(trap3, Point3{T}(50,-50,0), vy) ≈ 50
    @test distanceToIn(trap3, Point3{T}(50,-50,0), vmy) ≈ Inf
    @test distanceToIn(trap4, Point3{T}(50, 50,0), vy) ≈ Inf
    @test distanceToIn(trap4, Point3{T}(50, 50,0), vmy) ≈ 50

    @test distanceToIn(trap1, Point3{T}(0,60,0), vxmy) ≈ Inf
    @test distanceToIn(trap1, Point3{T}(0,50,0), vxmy) ≈ √800
    @test distanceToIn(trap1, Point3{T}(0,40,0), vxmy) ≈ 10√2
    @test distanceToIn(trap1, Point3{T}(0,40,50), vxmy) ≈ Inf
  
    # Parallel to side planes

    @test distanceToIn(trap1, Point3{T}(40,60,0), vmx) ≈ Inf
    @test distanceToIn(trap1, Point3{T}(40,60,0), vmy) ≈ Inf
    @test distanceToIn(trap1, Point3{T}(40,60,50), vmz) ≈ Inf
    @test distanceToIn(trap1, Point3{T}(0,0,50), vymz) ≈ 10√2
    @test distanceToIn(trap1, Point3{T}(0,0,80), vymz) ≈ Inf
    @test distanceToIn(trap1, Point3{T}(0,0,70), vymz) ≈ 30√2

  
    # Check Extent and cached BBox
    low, high = extent(trap1) 
    @test low  ≈ Point3{T}(-20,-30,-40)
    @test high ≈ Point3{T}( 20, 30, 40)

    low, high = extent(trap2)
    @test low  ≈ Point3{T}(-30,-40,-40)
    @test high ≈ Point3{T}( 30, 40, 40)

end


