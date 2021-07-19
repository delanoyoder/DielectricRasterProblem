s = Square(0.4, 1.1, 0.25, 1.0)

u = SVector(0.0, 1.0); d = -u
r = SVector(1.0, 0.0); l = -r
ul = SVector(-1/sqrt(2), 1/sqrt(2)); dr = -ul
ur = SVector(1/sqrt(2), 1/sqrt(2)); dl = -ur

@testset "External" begin
    @test find_nearest_surface_normal(s, (0,0)) ≈ dl
    @test find_nearest_surface_normal(s, (0,1.0)) ≈ l
    @test find_nearest_surface_normal(s, (0,1.23)) ≈ ul
    @test find_nearest_surface_normal(s, (0.523,1.99)) ≈ u
    @test find_nearest_surface_normal(s, (0.526,1.226)) ≈ ur
    @test find_nearest_surface_normal(s, (0.526,1.1)) ≈ r
    @test find_nearest_surface_normal(s, (0.526,0.974)) ≈ dr
    @test find_nearest_surface_normal(s, (0.4,0.9)) ≈ d
end

@testset "Internal" begin
    @test find_nearest_surface_normal(s, (0.35,1.05)) ≈ dl
    @test find_nearest_surface_normal(s, (0.35,1.1)) ≈ l
    @test find_nearest_surface_normal(s, (0.35,1.15)) ≈ ul
    @test find_nearest_surface_normal(s, (0.4,1.15)) ≈ u
    @test find_nearest_surface_normal(s, (0.45,1.15)) ≈ ur
    @test find_nearest_surface_normal(s, (0.45,1.1)) ≈ r
    @test find_nearest_surface_normal(s, (0.45,1.05)) ≈ dr
    @test find_nearest_surface_normal(s, (0.4,1.05)) ≈ d
end
