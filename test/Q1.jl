@testset "Rastering Square" begin
    A = GriddedArray((0.0, 1.0), (0.0, 2.0), fill(1.0, 5, 5))
    s = Square(0.5, 1.0, 0.75, 1.0)
    pixel_area = (1/4) * (2/4)

    raster_area!(A, s)

    # Test consistency of raster
    @test sum(A.pixels * pixel_area) == (s.w)^2

    # Test correctness of raster
    raster = zeros(5,5)
    raster[2:4,[2, 4]] .= 0.25; raster[2:4, 3] .= 1.0

    @test A.pixels ≈ raster
end
