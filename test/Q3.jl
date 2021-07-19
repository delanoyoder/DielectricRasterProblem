@testset "Raster Scalar" begin
    # Identity
    A = GriddedArray((0.0, 1.0), (0.0, 2.0), fill(1.0, 5, 5))
    s = Square(0.5, 1.0, 0.75, 1.0)

    raster_harmonic!(A, s)
    @test A.pixels ≈ ones(5,5)

    # Varying
    A = GriddedArray((0.0, 1.0), (0.0, 2.0), fill(2.0, 5, 5))
    s = Square(0.5, 1.0, 0.75, 1.5)
    result = fill(2.0, 5, 5)
    result[2:4, [2,4]] .= 24/13
    result[2:4, 3] .= 1.5

    raster_harmonic!(A, s)
    @test A.pixels ≈ result
end
