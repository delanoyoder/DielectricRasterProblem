@testset "Raster Tensor" begin
    output = TensorGriddedArray((0.0, 2.0), (0.0, 2.0), fill(SMatrix{2,2}(0.0, 0.0, 0.0, 0.0), 5, 5))
    dielectric = GriddedArray((0.0, 2.0), (0.0, 2.0), fill(2.0, 5, 5))
    square = Square(1.0, 1.0, 1.0, 1.5)

    raster_tensor_dielectric!(output, dielectric, square)

    # Testing SMatrices
    diagnl = SMatrix{2,2}(2.0, 0.0, 0.0, 2.0)
    ctr = SMatrix{2,2}(1.5, 0.0, 0.0, 1.5)
    side = SMatrix{2,2}(7/4, 0.0, 0.0, 12/7)

    a = 1.3258252147247767; b = 1.305427903729011
    crnr = SMatrix{2,2}(-a, -a, b, -b)

    R = SMatrix{2,2}(0.0, -1.0, 1.0, 0.0)

    # Border tensors
    for idxs in ((1, :), (5, :), (:, 1), (:, 5))
        @test all(output.pixels[idxs...] .== Ref(diagnl))
    end
    # Center tensor
    @test output.pixels[3,3] ≈ ctr
    # Sides
    @test output.pixels[3,4] ≈ side
    @test output.pixels[2,3] ≈ R*side
    @test output.pixels[3,2] ≈ R^2*side
    @test output.pixels[4,3] ≈ R^3*side
    # Corners
    @test output.pixels[2,2] ≈ crnr
    @test output.pixels[4,2] ≈ R*crnr
    @test output.pixels[4,4] ≈ R^2*crnr
    @test output.pixels[2,4] ≈ R^3*crnr
end
