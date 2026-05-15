using NAOMi
using Test

@testset "zernike_polynomial — low-order analytic forms" begin
    xs = -1.0:0.05:1.0
    X = [x for x in xs, _ in xs]
    Y = [y for _ in xs, y in xs]

    # Z_1 = piston (constant 1)
    Z1 = zernike_polynomial(1, X, Y)
    @test all(Z1 .≈ 1.0)

    # Z_2 = 2 r cos(θ) = 2 x  (NB: upstream's iseven(j) branch picks cos)
    Z2 = zernike_polynomial(2, X, Y)
    @test maximum(abs.(Z2 .- 2 .* X)) < 1e-12

    # Z_3 = 2 r sin(θ) = 2 y
    Z3 = zernike_polynomial(3, X, Y)
    @test maximum(abs.(Z3 .- 2 .* Y)) < 1e-12

    # Z_4 = sqrt(3) (2 r^2 - 1) — defocus
    Z4 = zernike_polynomial(4, X, Y)
    expected = sqrt(3) .* (2 .* (X .^ 2 .+ Y .^ 2) .- 1)
    @test maximum(abs.(Z4 .- expected)) < 1e-12
end

@testset "zernike_polynomial — orthonormality on unit disk" begin
    h = 0.01
    xs = -1.0:h:1.0
    X = [x for x in xs, _ in xs]
    Y = [y for _ in xs, y in xs]
    disk = X .^ 2 .+ Y .^ 2 .≤ 1.0
    # Each diagonal entry ≈ π; off-diagonals ≈ 0 (numerical integration on
    # a Cartesian grid converges modestly — keep tolerances loose).
    for j in (1, 4, 11)
        Zj = zernike_polynomial(j, X, Y)
        integ = sum(Zj .^ 2 .* disk) * h^2
        @test isapprox(integ, π; rtol = 0.05)
    end
    for (i, j) in ((1, 4), (2, 3), (4, 11), (5, 6))
        ij = sum(zernike_polynomial(i, X, Y) .* zernike_polynomial(j, X, Y) .* disk) * h^2
        @test abs(ij) < 0.05
    end
end

@testset "zernike_polynomial — argument validation" begin
    @test_throws ArgumentError zernike_polynomial(1, zeros(3, 3), zeros(2, 2))
    @test_throws ArgumentError zernike_polynomial(0, zeros(2, 2), zeros(2, 2))
end

@testset "generate_zernike_weights — λ scaling" begin
    pp = PSFParams()
    abb = generate_zernike_weights(pp)
    @test length(abb) == length(pp.zernikeWt)
    @test all(abb .≈ pp.zernikeWt .* pp.lambda .* 1e-6)
end

@testset "apply_zernike — phase-only preserves |U|²" begin
    xs = -1.0:0.05:1.0
    X = [x for x in xs, _ in xs]
    Y = [y for _ in xs, y in xs]
    Uin = ones(ComplexF64, size(X))
    Uout = apply_zernike(Uin, X, Y, 1.0, [0.0, 0.0, 0.0, 0.1])
    @test maximum(abs.(abs2.(Uout) .- 1.0)) < 1e-12

    # Zero aberration → identity (up to a global mean-subtraction)
    U0 = apply_zernike(Uin, X, Y, 1.0, zeros(5))
    @test maximum(abs.(U0 .- Uin)) < 1e-12

    # Shape mismatch
    @test_throws ArgumentError apply_zernike(zeros(ComplexF64, 4, 4),
                                              zeros(3, 3), zeros(3, 3),
                                              1.0, [0.1])
end

@testset "generate_back_aperture — shape and aperture mask" begin
    vp = VolumeParams(vol_sz = [20, 20, 10], vres = 2.0)
    finalize!(vp)
    pp = PSFParams()
    U = generate_back_aperture(vp, pp)
    @test size(U, 1) == size(U, 2)
    @test size(U, 1) > 0
    @test all(isfinite, U)
    # Outside the hard aperture (set by objrad), |U| is zero.
    n = size(U, 1)
    @test U[1, 1] == 0
    @test abs(U[n ÷ 2, n ÷ 2]) > 0
end

@testset "generate_back_aperture — zero aberrations → real-aperture profile" begin
    vp = VolumeParams(vol_sz = [20, 20, 10], vres = 2.0)
    finalize!(vp)
    pp = PSFParams()
    pp.zernikeWt = zeros(length(pp.zernikeWt))
    U = generate_back_aperture(vp, pp)
    # With zero abb, apply_zernike is skipped; the field should equal the
    # raw Gaussian profile (real, after multiplication by the focusing
    # phase). We can't easily reconstruct it without re-running the
    # internals, so just sanity-check the structure: |U| peaks at the
    # centre and falls off radially.
    n = size(U, 1)
    c = n ÷ 2
    @test abs(U[c, c]) ≥ abs(U[c + 50, c])
end
