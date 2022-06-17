using Test
using Geom4hep, LinearAlgebra

include("../src/Units.jl")

@testset "Polycone{$T}" for T in (Float64, Float32)
    # Polycones
    myPCone = Polycone{T}([30, 30, 0, 0, 0, 0, 40, 40],[70, 70, 70, 40, 40, 80, 80, 60],[-20, -10, -10, 0, 10, 20, 30, 40], -10deg, 10deg)
    simple = Polycone{T}([0,0,0], [70,70,80], [-10,0,10], 0, 2π)
    cms_trak = Polycone{T}([74., 34., 31., 31., 31., 31., 34., 74.],
                           [1233., 1233., 1233., 1233., 1233., 1233., 1233., 1233.],
                           [-2935., -1899., -1899., -1899., 1899., 1899., 1899., 2935.], 0., 2π)

    Nz = 4
    rmin = T[0.1, 0.0, 0.0, 0.4]
    rmax = T[1., 2., 2., 1.5]
    z    = T[-1, -0.5, 0.5, 2]
  
    poly1 = Polycone{T}(rmin, rmax, z, 0., 2π)
    # let's make external separate cones representing the section
    section1 = Cone{T}(rmin[1], rmax[1], rmin[2], rmax[2], (z[2] - z[1]) / 2., 0, 2π)
    section2 = Cone{T}(rmin[2], rmax[2], rmin[3], rmax[3], (z[3] - z[2]) / 2., 0, 2π)
    section3 = Cone{T}(rmin[3], rmax[3], rmin[4], rmax[4], (z[4] - z[3]) / 2., 0, 2π)

    @test getNz(poly1) == 4
    @test getNSections(poly1) == 3
    @test getSectionIndex(poly1, -0.8) == 1
    @test getSectionIndex(poly1, 0.51) == 3
    @test getSectionIndex(poly1, 0.) == 2
    @test getSectionIndex(poly1, -2.) == -1
    @test getSectionIndex(poly1, 3.) == -2

    @test capacity(poly1) > 0.
    @test capacity(poly1) ≈ capacity(section1) + capacity(section2) + capacity(section3)

    # test contains/inside
    @test contains(poly1, Point3{T}(0., 0., 0.)) == true
    @test contains(poly1, Point3{T}(0., 0., -2.)) == false
    @test contains(poly1, Point3{T}(0., 0., -0.8)) == false
    @test contains(poly1, Point3{T}(0., 0., -1.8)) == false
    @test contains(poly1, Point3{T}(0., 0., 10.)) == false
    @test contains(poly1, Point3{T}(0., 0., 1.8)) == false

    # test DistanceToIn
    @test distanceToIn(poly1, Point3{T}(0, 0, -3), Vector3{T}(0, 0, 1)) == 2.5
    @test distanceToIn(poly1, Point3{T}(0, 0, -3), Vector3{T}(0, 0, -1)) == Inf
    @test distanceToIn(poly1, Point3{T}(0, 0, 3), Vector3{T}(0, 0, -1)) == 2.5
    @test distanceToIn(poly1, Point3{T}(0, 0, 3), Vector3{T}(0, 0, 1)) == Inf
    @test distanceToIn(poly1, Point3{T}(3, 0, 0), Vector3{T}(-1, 0, 0)) == 1.
    @test isapprox(distanceToIn(poly1, Point3{T}(0, 0, 2 - 10*kTolerance(T)), Vector3{T}(1, 0, 0)), 0.4, atol=10*kTolerance(T))

    # test SafetyToIn
    @test safetyToIn(poly1, Point3{T}(0, 0, -3)) == 2.
    @test safetyToIn(poly1, Point3{T}(0.5, 0, -1)) == 0.
    @test safetyToIn(poly1, Point3{T}(0, 0, 3)) == 1.
    @test safetyToIn(poly1, Point3{T}(2, 0, 0.1)) == 0.

    # test SafetyToOut
    @test safetyToOut(poly1, Point3{T}(0, 0, 0)) == 0.5
    @test safetyToOut(poly1, Point3{T}(0, 0, 0.5)) == 0.
    @test isapprox(safetyToOut(poly1, Point3{T}(1.9, 0, 0)), 0.1, atol=10*kTolerance(T))
    @test safetyToOut(poly1, Point3{T}(0.2, 0, -1)) == 0.
    @test safetyToOut(poly1, Point3{T}(1.4, 0, 2)) == 0.

    # test DistanceToOut
    @test distanceToOut(poly1, Point3{T}(0, 0, 0), Vector3{T}(0, 0, 1)) == 0.5
    @test distanceToOut(poly1, Point3{T}(0, 0, 0), Vector3{T}(0, 0, -1)) == 0.5
    @test distanceToOut(poly1, Point3{T}(2, 0, 0), Vector3{T}(1, 0, 0)) == 0.
    @test distanceToOut(poly1, Point3{T}(2, 0, 0), Vector3{T}(-1, 0, 0)) == 4.

    @test distanceToOut(poly1, Point3{T}(1, 0, 2), Vector3{T}(0, 0, 1)) == 0.
    @test distanceToOut(poly1, Point3{T}(0.5, 0, -1), Vector3{T}(0, 0, -1)) == 0.
    @test distanceToOut(poly1, Point3{T}(0.5, 0, -1), Vector3{T}(0, 0, 1)) == 3.

    # Check Cubic volume
    @test capacity(simple) ≈ T(π) * (70 * 70 * 10 + 10 * (70 * 70 + 80 * 80 + 70 * 80) / 3.)
    # Check surface area 
    @test surface(simple) ≈  T(π) * (70 * 70 + 80 * 80 + (70 + 80) * √(10 * 10 + 10 * 10) + 10 * 2 * 70)
    # Check Inside

    @test inside(simple, Point3{T}(0,0,0)) == kInside
    @test inside(simple, Point3{T}(0,0,0)) == kInside

    pzero = Point3{T}(0, 0, 0)
    ponxside = Point3{T}(70, 0, -5)
    ponyside = Point3{T}(0, 70, -5)
    ponzside = Point3{T}(70, 0, 10)
    ponmxside = Point3{T}(-70, 0, -5)
    ponmyside = Point3{T}(0, -70, -5)
    ponmzside = Point3{T}(0, 0, -10)
    ponzsidey = Point3{T}(0, 25, 0)
    ponmzsidey = Point3{T}(4, 25, 0)

    pbigx = Point3{T}(100, 0, 0)
    pbigy = Point3{T}(0, 100, 0)
    pbigz = Point3{T}(0, 0, 100)
    pbigmx = Point3{T}(-100, 0, 0)
    pbigmy = Point3{T}(0, -100, 0)
    pbigmz = Point3{T}(0, 0, -100)

    vx = Vector3{T}(1, 0, 0)
    vy = Vector3{T}(0, 1, 0)
    vz = Vector3{T}(0, 0, 1)
    vmx = Vector3{T}(-1, 0, 0)
    vmy = Vector3{T}(0, -1, 0)
    vmz = Vector3{T}(0, 0, -1)
    vxy = Vector3{T}(1/√2, 1/√2, 0)
    vmxy = Vector3{T}(-1/√2, 1/√2, 0)
    vmxmy = Vector3{T}(-1/√2, -1/√2, 0)
    vxmy = Vector3{T}(1/√2, -1/√2, 0)

    @test inside(simple, pzero) == kInside
    @test inside(simple, pbigz) == kOutside
    @test inside(simple, pbigx) == kOutside
    @test inside(simple, pbigy) == kOutside

    @test inside(simple, ponxside) == kSurface
    @test inside(simple, ponyside) == kSurface
    @test inside(simple, ponzside) == kSurface
    @test inside(simple, ponmxside) == kSurface
    @test inside(simple, ponmyside) == kSurface
    @test inside(simple, ponmzside) == kSurface
    @test inside(simple, ponmzsidey) == kInside

    # SafetyToOut(P)
    @test safetyToOut(simple, Point3{T}(5, 5, -5)) ≈ 5.
    @test safetyToOut(simple, Point3{T}(5, 5, 7)) ≈ 3.
    @test safetyToOut(simple, Point3{T}(69, 0, -5)) ≈ 1.
    @test safetyToOut(simple, Point3{T}(-3, -3, 8)) ≈ 2.

    # DistanceToOut(P,V)
    @test distanceToOut(simple, pzero, vx) ≈ 70.
    @test distanceToOut(simple, pzero, vmx) ≈ 70.
    @test distanceToOut(simple, pzero, vy) ≈ 70.
    @test distanceToOut(simple, pzero, vmy) ≈ 70.
    @test distanceToOut(simple, pzero, vz) ≈ 10.
    @test distanceToOut(simple, Point3{T}(70, 0, -10), vx) ≈ 0.
    @test distanceToOut(simple, Point3{T}(-70, 0, -1), vmx) ≈ 0.
    @test distanceToOut(simple, Point3{T}(0, 70, -10), vy) ≈ 0.
    @test distanceToOut(simple, Point3{T}(0, -70, -1), vmy) ≈ 0.

    # SafetyToIn(P)
    @test safetyToIn(simple, pbigx) ≈ T(21.213203435596423)
    @test safetyToIn(simple, pbigmx) ≈ T(21.213203435596423)
    @test safetyToIn(simple, pbigy) ≈ T(21.213203435596423)
    @test safetyToIn(simple, pbigmy) ≈ T(21.213203435596423)
    @test safetyToIn(simple, pbigz) ≈ 90.
    @test safetyToIn(simple, pbigmz) ≈ 90.

    # DistanceToIn(P,V)
    @test distanceToIn(simple, Point3{T}(100, 0, -1), vmx) ≈ 30.
    @test distanceToIn(simple, Point3{T}(-100, 0, -1), vx) ≈ 30.
    @test distanceToIn(simple, Point3{T}(0, 100, -5), vmy) ≈ 30.
    @test distanceToIn(simple, Point3{T}(0, -100, -5), vy) ≈ 30.
    @test distanceToIn(simple, pbigz, vmz) ≈ 90.
    @test distanceToIn(simple, pbigmz, vz) ≈ 90.
    @test distanceToIn(simple, pbigx, vxy) == Inf
    @test distanceToIn(simple, pbigmx, vmxy) == Inf

    # Extent
    low, high = extent(simple)
    @test low  ≈ Point3{T}(-80.0, -80.0, -10.0)
    @test high ≈ Point3{T}(80.0, 80.0, 10.0)
end  

#=

  // check that Normal() returns valid=false and a non-zero normal for points away from the surface

  Vec_t point(70, 70, -5);
  if ((valid = Simple.Normal(point, normal)) || !ApproxEqual<Precision>(normal.Mag2(), 1))
    std::cout << "Simple.Normal() normal not checked: Line " << __LINE__ << ", p=" << point << ", normal=" << normal
              << ", valid=" << valid << "\n";
  point.z() = -10;
  if ((valid = Simple.Normal(point, normal)) || !ApproxEqual<Precision>(normal.Mag2(), 1))
    std::cout << "Simple.Normal() normal not checked: Line " << __LINE__ << ", p=" << point << ", normal=" << normal
              << ", valid=" << valid << "\n";
  if ((valid = Simple.Normal(pbigz, normal)) || !ApproxEqual<Precision>(normal.Mag2(), 1))
    std::cout << "Simple.Normal() normal not checked: Line " << __LINE__ << ", p=" << pbigz << ", normal=" << normal
              << ", valid=" << valid << "\n";
  if ((valid = Simple.Normal(pbigmz, normal)) || !ApproxEqual<Precision>(normal.Mag2(), 1))
    std::cout << "Simple.Normal() normal not checked: Line " << __LINE__ << ", p=" << pbigmz << ", normal=" << normal
              << ", valid=" << valid << "\n";

  // Check Surface Normal

  valid = Simple.Normal(ponxside, normal);
  assert(ApproxEqual(normal, Vec_t(1, 0, 0)) && valid);
  valid = Simple.Normal(ponmxside, normal);
  assert(ApproxEqual(normal, Vec_t(-1, 0, 0)));
  valid = Simple.Normal(ponyside, normal);
  assert(ApproxEqual(normal, Vec_t(0, 1, 0)));
  valid = Simple.Normal(Vec_t(0, 0, 10), normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));
  valid = Simple.Normal(Vec_t(0, 0, -10), normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));

  // Normals on Edges

  Vec_t edgeXZ(80.0, 0.0, 10.0);
  Vec_t edgeYZ(0., 80.0, 10.0);
  Vec_t edgeXmZ(70.0, 0.0, -10.0);
  Vec_t edgeYmZ(0.0, 70.0, -10.0);
  Vec_t edgemXZ(-80.0, 0.0, 10.0);
  Vec_t edgemYZ(0., -80.0, 10.0);
  Vec_t edgemXmZ(-70.0, 0.0, -10.0);
  Vec_t edgemYmZ(0.0, -70.0, -10.0);
  // Precision invSqrt2 = 1.0 / std::sqrt(2.0);
  // Precision invSqrt3 = 1.0 / std::sqrt( 3.0);

  valid = Simple.Normal(edgeXmZ, normal);
  // assert(ApproxEqual(normal, Vec_t(invSqrt2, 0.0, -invSqrt2)));
  valid = Simple.Normal(edgemXmZ, normal);
  // assert(ApproxEqual(normal, Vec_t(-invSqrt2, 0.0, -invSqrt2)));
  valid = Simple.Normal(edgeYmZ, normal);
  // assert(ApproxEqual(normal, Vec_t(0.0, invSqrt2, -invSqrt2)));
  valid = Simple.Normal(edgemYmZ, normal);
  // assert(ApproxEqual(normal, Vec_t(0.0, -invSqrt2, -invSqrt2)));

  const Precision xyn = 0.92388, zn = 0.382683;
  valid = Simple.Normal(edgeXZ, normal);
  std::cout << "Simple.Normal(): p=" << edgeXZ << ", normal=" << normal << ", valid=" << valid << std::endl;
  assert(ApproxEqual(normal, Vec_t(xyn, 0, zn)));
  valid = Simple.Normal(edgemXZ, normal);
  assert(ApproxEqual(normal, Vec_t(-xyn, 0, zn)));
  valid = Simple.Normal(edgeYZ, normal);
  assert(ApproxEqual(normal, Vec_t(0, xyn, zn)));
  valid = Simple.Normal(edgemYZ, normal);
  assert(ApproxEqual(normal, Vec_t(0, -xyn, zn)));

 
  // Check Extent and cached BBox
  Vec_t minExtent, maxExtent;
  Vec_t minBBox, maxBBox;
  Simple.Extent(minExtent, maxExtent);
  Simple.GetUnplacedVolume()->GetBBox(minBBox, maxBBox);
  // std::cout<<" min="<<minExtent<<" max="<<maxExtent<<std::endl;
  assert(ApproxEqual(minExtent, Vec_t(-80, -80, -10)));
  assert(ApproxEqual(maxExtent, Vec_t(80, 80, 10)));
  assert(ApproxEqual(minExtent, minBBox));
  assert(ApproxEqual(maxExtent, maxBBox));
  MyPCone->Extent(minExtent, maxExtent);
  MyPCone->GetUnplacedVolume()->GetBBox(minBBox, maxBBox);
  // std::cout<<" min="<<minExtent<<" max="<<maxExtent<<std::endl;
  // assert(ApproxEqual(minExtent, Vec_t(-80, -80, -20)));
  // assert(ApproxEqual(maxExtent, Vec_t(80, 80, 40)));
  assert(ApproxEqual(minExtent, minBBox));
  assert(ApproxEqual(maxExtent, maxBBox));
=#