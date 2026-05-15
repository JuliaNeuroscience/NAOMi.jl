using NAOMi
using Test

@testset "gaussian_psf — shape and central value" begin
    psflen, lambda, sampling, matSize = 4.0, 0.92, [0.1, 0.1, 0.1], [51, 51, 81]
    psf, x, y, z = gaussian_psf(psflen, lambda, sampling, matSize)
    @test size(psf) == (51, 51, 81)
    @test length(x) == 51 && length(y) == 51 && length(z) == 81
    # Coordinate axes are centered on zero (using round(N/2))
    @test x[round(Int, 51 / 2)] == 0.0
    @test y[round(Int, 51 / 2)] == 0.0
    @test z[round(Int, 81 / 2)] == 0.0
    # PSF is real, positive, max at origin
    @test all(psf .≥ 0)
    @test psf[round(Int, 51 / 2), round(Int, 51 / 2), round(Int, 81 / 2)] ≈ 1.0
    @test maximum(psf) ≈ 1.0
end

@testset "gaussian_psf_na — NA scaling and half-power point" begin
    na, lambda, nidx = 0.6, 0.92, 1.33
    sampling, matSize = [0.1, 0.1, 0.1], [51, 51, 401]
    psflen = 0.626 * lambda / (nidx - sqrt(nidx^2 - na^2))
    psf, x, y, z, intensity = gaussian_psf_na(na, lambda, sampling, matSize;
                                              nidx = nidx)
    # Center indexing for this function uses 0-based offset + round(N/2);
    # the center is at round(Int, N/2) + 1.
    cx = round(Int, 51 / 2) + 1
    cz = round(Int, 401 / 2) + 1
    @test intensity[cx, cx, cz] ≈ 1.0
    # Axial half-power: intensity at z = psflen/2 should be ≈ 0.5
    z_half_offset = round(Int, (psflen / 2) / sampling[3])
    @test isapprox(intensity[cx, cx, cz + z_half_offset], 0.5; atol = 0.01)
    @test isapprox(intensity[cx, cx, cz - z_half_offset], 0.5; atol = 0.01)
    # psf is intensity squared
    @test isapprox(psf[cx, cx, cz + z_half_offset], 0.25; atol = 0.02)
    # Sweep NA: doubling NA shortens psflen
    p1 = 0.626 * lambda / (1.33 - sqrt(1.33^2 - 0.3^2))
    p2 = 0.626 * lambda / (1.33 - sqrt(1.33^2 - 0.6^2))
    @test p1 > p2
end

@testset "gaussian_beam_size — monotone in dist, zero axial" begin
    p = PSFParams()
    sz1 = gaussian_beam_size(p, 50.0)
    sz2 = gaussian_beam_size(p, 100.0)
    @test sz1[1] > 0 && sz1[2] > 0 && sz1[3] == 0
    @test sz1[1] == sz1[2]
    @test sz2[1] > sz1[1]
    # apod scales linearly
    sz3 = gaussian_beam_size(p, 100.0; apod = 4)
    @test sz3[1] == 2 * sz2[1]
end

@testset "generate_gaussian_profile — aperture and centring" begin
    n = 41
    rng = -1.0:0.05:1.0                    # n = 41 samples
    X = [x for x in rng, _ in rng]
    Y = [y for _ in rng, y in rng]
    k = 2π / 0.92
    rad, aper, fl = 0.2, 0.5, 4.5
    U = generate_gaussian_profile(X, Y, rad, aper, k, fl)
    @test size(U) == (n, n)
    @test all(isfinite, U)
    # Outside the aperture the field is zero
    outside = @. (X^2 + Y^2 ≥ aper^2)
    @test all(U[outside] .== 0)
    # Maximum |U| at the Gaussian centre (offset = (0, 0))
    @test argmax(abs.(U)) == CartesianIndex((n + 1) ÷ 2, (n + 1) ÷ 2)
    # With nonzero offset, max moves
    off = (0.3, 0.0)
    U2 = generate_gaussian_profile(X, Y, rad, aper, k, fl, off)
    idx = argmax(abs.(U2))
    @test X[idx] > 0
end

@testset "generate_gaussian_profile — shape validation" begin
    X = zeros(5, 5); Y = zeros(5, 6)
    @test_throws ArgumentError generate_gaussian_profile(X, Y, 0.1, 0.5, 1.0, 1.0)
end
