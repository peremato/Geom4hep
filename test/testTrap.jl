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
    @test safety_to_out(trap1, pzero) ≈ 20
    @test safety_to_out(trap1, Point3{T}(1,0,0)) ≈ 19
    @test safety_to_out(trap1, Point3{T}(0,1,0)) ≈ 20
    @test safety_to_out(trap1, Point3{T}(0,0,1)) ≈ 20

    @test safety_to_out(trap2, pzero) ≈ 20*4/√17
    @test safety_to_out(trap2, Point3{T}(1,0,0)) ≈ 19*4/√17
    @test safety_to_out(trap2, Point3{T}(0,1,0)) ≈ 20*4/√17
    @test safety_to_out(trap2, Point3{T}(0,0,1)) ≈ 20*4/√17 + 1/√17

    # SafetyToIn(P)
    @test safety_to_in(trap1, pbig) ≈ 80
    @test safety_to_in(trap1, pbigx) ≈ 80
    @test safety_to_in(trap1, pbigmx) ≈ 80
    @test safety_to_in(trap1, pbigy) ≈ 70
    @test safety_to_in(trap1, pbigmy) ≈ 70
    @test safety_to_in(trap1, pbigz) ≈ 60
    @test safety_to_in(trap1, pbigmz) ≈ 60

    @test safety_to_in(trap2, pbigx) ≈ 80*4/√17
    @test safety_to_in(trap2, pbigmx) ≈ 80*4/√17
    @test safety_to_in(trap2, pbigy) ≈ 70*4/√17
    @test safety_to_in(trap2, pbigmy) ≈ 70*4/√17
    @test safety_to_in(trap2, pbigz) ≈ 60
    @test safety_to_in(trap2, pbigmz) ≈ 60


    # DistanceToOut(P,V)
    @test distance_to_out(trap1, pzero, vx) ≈ 20
    @test distance_to_out(trap1, pzero, vmx) ≈ 20
    @test distance_to_out(trap1, pzero, vy) ≈ 30
    @test distance_to_out(trap1, pzero, vmy) ≈ 30
    @test distance_to_out(trap1, pzero, vz) ≈ 40
    @test distance_to_out(trap1, pzero, vmz) ≈ 40
    @test distance_to_out(trap1, pzero, vxy) ≈ √800
    @test distance_to_out(trap1, ponxside, vx) ≈ 0
    @test distance_to_out(trap1, ponmxside, vmx) ≈ 0
    @test distance_to_out(trap1, ponyside, vy) ≈ 0
    @test distance_to_out(trap1, ponmyside, vmy) ≈ 0
    @test distance_to_out(trap1, ponzside, vz) ≈ 0
    @test distance_to_out(trap1, ponmzside, vmz) ≈ 0
    @test distance_to_out(trap1, ponxside, vmx) ≈ 40
    @test distance_to_out(trap1, ponmxside, vx) ≈ 40
    @test distance_to_out(trap1, ponyside, vmy) ≈ 60
    @test distance_to_out(trap1, ponmyside, vy) ≈ 60
    @test distance_to_out(trap1, ponzside, vmz) ≈ 80
    @test distance_to_out(trap1, ponmzside, vz) ≈ 80

    @test distance_to_out(trap2, pzero, vx) ≈ 20
    @test distance_to_out(trap2, pzero, vmx) ≈ 20
    @test distance_to_out(trap2, pzero, vy) ≈ 30
    @test distance_to_out(trap2, pzero, vmy) ≈ 30
    @test distance_to_out(trap2, pzero, vz) ≈ 40
    @test distance_to_out(trap2, pzero, vmz) ≈ 40
    @test distance_to_out(trap2, pzero, vxy) ≈ √800
    @test distance_to_out(trap2, ponxside, vx) ≈ 0
    @test distance_to_out(trap2, ponmxside, vmx) ≈ 0
    @test isapprox(distance_to_out(trap2, ponyside, vy), 0, atol=1e-6)
    @test distance_to_out(trap2, ponmyside, vmy) ≈ 0
    @test distance_to_out(trap2, ponzside, vz) ≈ 0
    @test distance_to_out(trap2, ponmzside, vmz) ≈ 0

    # DistanceToIn(P,V)
    @test distance_to_in(trap1, pbigx, vmx) ≈ 80
    @test distance_to_in(trap1, pbigmx, vx) ≈ 80
    @test distance_to_in(trap1, pbigy, vmy) ≈ 70
    @test distance_to_in(trap1, pbigmy, vy) ≈ 70
    @test distance_to_in(trap1, pbigz, vmz) ≈ 60
    @test distance_to_in(trap1, pbigmz, vz) ≈ 60
    @test distance_to_in(trap1, pbigx, vxy) ≈ Inf
    @test distance_to_in(trap1, pbigmx, vxy) ≈ Inf

    @test distance_to_in(trap2, pbigx, vmx) ≈ 80
    @test distance_to_in(trap2, pbigmx, vx) ≈ 80
    @test distance_to_in(trap2, pbigy, vmy) ≈ 70
    @test distance_to_in(trap2, pbigmy, vy) ≈ 70
    @test distance_to_in(trap2, pbigz, vmz) ≈ 60
    @test distance_to_in(trap2, pbigmz, vz) ≈ 60
    @test distance_to_in(trap2, pbigx, vxy) ≈ Inf
    @test distance_to_in(trap2, pbigmx, vxy) ≈ Inf

    @test distance_to_in(trap3, Point3{T}(50,-50,0), vy) ≈ 50
    @test distance_to_in(trap3, Point3{T}(50,-50,0), vmy) ≈ Inf
    @test distance_to_in(trap4, Point3{T}(50, 50,0), vy) ≈ Inf
    @test distance_to_in(trap4, Point3{T}(50, 50,0), vmy) ≈ 50

    @test distance_to_in(trap1, Point3{T}(0,60,0), vxmy) ≈ Inf
    @test distance_to_in(trap1, Point3{T}(0,50,0), vxmy) ≈ √800
    @test distance_to_in(trap1, Point3{T}(0,40,0), vxmy) ≈ 10√2
    @test distance_to_in(trap1, Point3{T}(0,40,50), vxmy) ≈ Inf
  
    # Parallel to side planes

    @test distance_to_in(trap1, Point3{T}(40,60,0), vmx) ≈ Inf
    @test distance_to_in(trap1, Point3{T}(40,60,0), vmy) ≈ Inf
    @test distance_to_in(trap1, Point3{T}(40,60,50), vmz) ≈ Inf
    @test distance_to_in(trap1, Point3{T}(0,0,50), vymz) ≈ 10√2
    @test distance_to_in(trap1, Point3{T}(0,0,80), vymz) ≈ Inf
    @test distance_to_in(trap1, Point3{T}(0,0,70), vymz) ≈ 30√2

  
    # Check Extent and cached BBox
    low, high = extent(trap1) 
    @test low  ≈ Point3{T}(-20,-30,-40)
    @test high ≈ Point3{T}( 20, 30, 40)

    low, high = extent(trap2)
    @test low  ≈ Point3{T}(-30,-40,-40)
    @test high ≈ Point3{T}( 30, 40, 40)

end



@testset "Trap with degenerate faces{$T}" for T in (Float64, Float32)
    # a face at +z with sides of zero length (a triangular prism), as the TEC module frames of CMS
    tec = Trap{T}(67.15, 0.0494756, 0, 0.382, 6.65, 6.65, 0, 0.382, 5e-14, 5e-14, 0)
    wedge = Trap{T}(20, 0, 0, 10, 10, 10, 0, 10, 0, 0, 0)         # edge at +z along y
    @test all(p -> all(isfinite, p.normal) && norm(p.normal) ≈ 1, tec.planes)
    @test all(p -> all(isfinite, p.normal) && norm(p.normal) ≈ 1, wedge.planes)
    # entering through the face at -y, leaving through the one at +y (it used to return 1.2e9)
    p, d = Point3{T}(-4.451195, -0.38195857, -13.122291), Vector3{T}(0, 1, -4.371139f-8)
    @test inside(tec, p) != kOutside                               # 4e-5 inside the -y face
    @test distance_to_out(tec, p, d) ≈ 0.764 rtol = 1e-3
    @test distance_to_in(wedge, Point3{T}(0, -50, 0), Vector3{T}(0, 1, 0)) ≈ 40
    @test distance_to_out(wedge, Point3{T}(0, 0, 0), Vector3{T}(0, 1, 0)) ≈ 10
    @test distance_to_out(wedge, Point3{T}(0, 0, 0), Vector3{T}(1, 0, 0)) ≈ 5      # side going to the edge
    @test capacity(wedge) ≈ 20 * 20 * 40 / 2
    # random rays: end points on the surface, and the shape entered or left there
    rng = MersenneTwister(11)
    for s in (tec, wedge)
        lo, hi = extent(s)
        nbad = 0
        for _ in 1:2000
            q = Point3{T}(lo + (hi - lo) .* (rand(rng, T, 3) .* 1.4 .- 0.2))
            v = normalize(Vector3{T}(randn(rng, T, 3)))
            pos = inside(s, q)
            if pos == kInside
                t = distance_to_out(s, q, v)
                e = q + t * v
                nbad += !(isfinite(t) && inside(s, e) == kSurface && inside(s, e + 10tolerance(T) * v) != kInside)
            elseif pos == kOutside
                t = distance_to_in(s, q, v)
                t == Inf && continue
                e = q + t * v
                nbad += !(inside(s, e) == kSurface && inside(s, e + 10tolerance(T) * v) != kOutside)
            end
        end
        @test nbad <= 5
    end
end
