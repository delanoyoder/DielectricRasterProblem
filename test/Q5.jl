e1 = SVector(1.0, 0.0)
e2 = SVector(0.0, 1.0)

@testset "Vectors" begin
    @test orient_along_normal(e1, e2) ≈ e1
    @test orient_along_normal(e1, e1) ≈ e2
    @test orient_along_normal(e2, e1) ≈ -e1
    @test orient_along_normal(e1 + e2, e1) ≈ e2 - e1
    @test orient_along_normal(
        cos(π/3)*e1 + sin(π/3)*e2, (e1 + e2)/sqrt(2)) ≈ (cos(7π/12)*e1 + sin(7π/12)*e2
        )
    @test orient_along_normal(
        -8*e1 + 7*e2, cos(3π/4)*e1 + sin(3π/4)*e2
        ) ≈ -8*cos(π/4)*(e1 - e2) + 7*sin(π/4)*(e2 + e1)
end

@testset "Matrices" begin
    @test orient_along_normal(hcat(e1, e2), e2) ≈ hcat(e1, e2)
    @test orient_along_normal(hcat(e2, e1), e1) ≈ hcat(-e1, e2)
    @test orient_along_normal(
        hcat((e1 + e2)/2, (e1 - e2)/2), -(e1 + e2)/sqrt(2)
        ) ≈ hcat(-e2/sqrt(2), -e1/sqrt(2))
end
