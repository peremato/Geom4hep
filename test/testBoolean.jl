@testset "Boolean{$T}" for T in (Float64, Float32)
    b10 = Box{T}(10, 10, 10)
    b5  = Box{T}(5, 5, 5)
    shiftx = Transformation3D{T}(10, 0, 0)   # right operand centred at x = +10
    xdir  = Vector3{T}(1, 0, 0)
    mxdir = Vector3{T}(-1, 0, 0)
    P(x, y=0, z=0) = Point3{T}(x, y, z)

    @testset "Union" begin
        u = BooleanUnion(b10, b10, shiftx)     # x ∈ [-10, 20]
        lo, hi = extent(u)
        @test lo ≈ P(-10, -10, -10)
        @test hi ≈ P(20, 10, 10)
        @test inside(u, P(0)) == kInside
        @test inside(u, P(15)) == kInside
        @test inside(u, P(10)) == kInside      # inside both operands
        @test inside(u, P(20)) == kSurface
        @test inside(u, P(25)) == kOutside
        @test inside(u, P(0, 15)) == kOutside
        @test distance_to_in(u, P(-30), xdir) ≈ 20
        @test distance_to_in(u, P(40), mxdir) ≈ 20
        @test distance_to_out(u, P(0), xdir) ≈ 20 atol = 10push_tolerance(T)
        @test distance_to_out(u, P(15), mxdir) ≈ 25 atol = 10push_tolerance(T)
    end

    @testset "Subtraction" begin
        s = BooleanSubtraction(b10, b5)        # hollow box: 5 < |x|,|y|,|z| ... < 10
        lo, hi = extent(s)
        @test lo ≈ P(-10, -10, -10)
        @test hi ≈ P(10, 10, 10)
        @test inside(s, P(0)) == kOutside
        @test inside(s, P(7)) == kInside
        @test inside(s, P(5)) == kSurface
        @test inside(s, P(10)) == kSurface
        @test inside(s, P(12)) == kOutside
        @test distance_to_in(s, P(-30), xdir) ≈ 20
        @test distance_to_out(s, P(7), xdir) ≈ 3
        @test distance_to_out(s, P(7), mxdir) ≈ 2
    end

    @testset "Intersection" begin
        i = BooleanIntersection(b10, b10, shiftx)  # x ∈ [0, 10]
        lo, hi = extent(i)
        @test lo ≈ P(0, -10, -10)
        @test hi ≈ P(10, 10, 10)
        @test inside(i, P(5)) == kInside
        @test inside(i, P(-5)) == kOutside
        @test inside(i, P(15)) == kOutside
        @test inside(i, P(0)) == kSurface
        @test distance_to_in(i, P(-30), xdir) ≈ 30
        @test distance_to_out(i, P(5), xdir) ≈ 5
        @test distance_to_out(i, P(5), mxdir) ≈ 5
    end

    @testset "Nested" begin
        # (b10 ∪ b10@+10) − b5 : hole around the origin
        n = BooleanSubtraction(BooleanUnion(b10, b10, shiftx), b5)
        @test inside(n, P(0)) == kOutside
        @test inside(n, P(7)) == kInside
        @test inside(n, P(15)) == kInside
        @test inside(n, P(25)) == kOutside
        @test distance_to_in(n, P(-30), xdir) ≈ 20
    end
end
