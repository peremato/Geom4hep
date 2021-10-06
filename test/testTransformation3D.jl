using Rotations

@testset "Trans3D{$T}" for T in (Float64, Float32)
    point = Point3{T}(-1., 1., 2.)
    t0 = one(Transformation3D{T})
    @test t0 == one(Transformation3D{T})
    @test isone(t0)
    @test !hasrotation(t0)
    @test !hastranslation(t0)
    @test transform(t0, point) == point
    @test t0 * point == point

    t1 = Transformation3D{T}(-2,-2,-2)
    @test !isone(t1)
    @test hastranslation(t1)
    @test !hasrotation(t1) 
    @test t1 * point == Point3{T}(1, 3, 4)


    t2 = Transformation3D{T}(2, 2, 2)
    @test t2 * (t1 * point) == point
    @test (t2 * t1) * point == point

    t3 = Transformation3D{T}(0, 0, 0, 0 , 0, π)
    @test !isone(t3)
    @test !hastranslation(t3)
    @test hasrotation(t3)
    @test t3 * point ≈ Point3{T}(1,-1, 2)

    t4 = Transformation3D{T}(0, 0, 1., RotXYZ{T}(0,0,π/4))
    @test t4 * Point3{T}(1,0,0) ≈ Point3{T}(√2/2, √2/2, -1)
    
    @test (t4 * Point3{T}(1,2,3)) * t4 ≈ Point3{T}(1,2,3)
end
