using NAOMi
using Test

@testset "fresnel_propagation_multi — vacuum focus of Gaussian back aperture" begin
    # Schmidt's matched-sampling rule for a back-to-focal-plane propagation.
    N = 128
    λ_m, nidx = 0.92e-6, 1.33
    λ_med = λ_m / nidx
    focal = 4.5e-3
    D2 = 0.5e-6
    D1 = λ_med * focal / (D2 * N)

    xs = ((-N ÷ 2):(N ÷ 2 - 1)) .* D1
    X = [x for x in xs, _ in xs]
    Y = [y for _ in xs, y in xs]
    rad = D1 * N / 4
    k_med = 2π / λ_med
    Uin = exp.(-((X .^ 2 .+ Y .^ 2) ./ rad^2)) .*
          exp.(-1im * k_med / (2 * focal) .* (X .^ 2 .+ Y .^ 2))
    phi = ones(ComplexF64, N, N, 2)

    Uout = fresnel_propagation_multi(Uin, λ_m, [D1, D2], [0.0, focal], phi, nidx)
    intensity = abs2.(Uout)

    # Output shape preserved
    @test size(Uout) == (N, N)
    # Focused spot is centered on (N÷2+1, N÷2+1).
    @test argmax(intensity) == CartesianIndex(N ÷ 2 + 1, N ÷ 2 + 1)
    # Lossless propagation: ∫|U|² dx² is conserved.
    energy_in  = sum(abs2.(Uin)) * D1^2
    energy_out = sum(abs2.(Uout)) * D2^2
    @test isapprox(energy_in, energy_out; rtol = 1e-6)
    # The 2-photon PSF (|U|⁴) is more peaked than the linear intensity.
    psf2p = intensity .^ 2
    @test argmax(psf2p) == argmax(intensity)
end

@testset "fresnel_propagation_multi — saveall returns per-plane stack" begin
    N = 32
    λ_m, nidx = 0.92e-6, 1.33
    Uin = ones(ComplexF64, N, N)
    phi = ones(ComplexF64, N, N, 3)
    z   = [0.0, 0.5e-4, 1.0e-4]
    dx  = [1e-6, 1e-6, 1e-6]
    Uout, Uall = fresnel_propagation_multi(Uin, λ_m, dx, z, phi, nidx;
                                            saveall = true)
    @test size(Uall) == (N, N, 3)
    # Last plane equals Uout (multiplied by final phase factor that is 1
    # when sc == 1 → no rescaling).
    @test Uall[:, :, end] ≈ Uout
end

@testset "fresnel_propagation_multi — argument validation" begin
    N = 16
    Uin = zeros(ComplexF64, N, N)
    phi = ones(ComplexF64, N, N, 2)
    # length(z) < 2
    @test_throws ArgumentError fresnel_propagation_multi(
        Uin, 1e-6, [1e-6], [0.0], phi, 1.33)
    # length(dx) ≠ length(z)
    @test_throws ArgumentError fresnel_propagation_multi(
        Uin, 1e-6, [1e-6], [0.0, 1e-3], phi, 1.33)
    # phi lateral shape mismatch
    bad_phi = ones(ComplexF64, N + 1, N, 2)
    @test_throws ArgumentError fresnel_propagation_multi(
        Uin, 1e-6, [1e-6, 1e-6], [0.0, 1e-3], bad_phi, 1.33)
    # size(phi, 3) ∉ {n, n-1}
    bad_phi3 = ones(ComplexF64, N, N, 5)
    @test_throws ArgumentError fresnel_propagation_multi(
        Uin, 1e-6, [1e-6, 1e-6], [0.0, 1e-3], bad_phi3, 1.33)
end

@testset "group_z_project — divisible and remainder paths" begin
    A = reshape(collect(1.0:24), 2, 3, 4)
    # Sum of two pairs of slices: each output entry = sum over groupsize=2.
    out_sum = group_z_project(A, 2; type = :sum)
    @test size(out_sum) == (2, 3, 2)
    @test out_sum[1, 1, 1] == A[1, 1, 1] + A[1, 1, 2]
    @test out_sum[2, 3, 2] == A[2, 3, 3] + A[2, 3, 4]
    # Mean of full-z reduces to a single plane.
    out_mean = group_z_project(A, 4; type = :mean)
    @test size(out_mean) == (2, 3, 1)
    @test out_mean[1, 1, 1] ≈ sum(A[1, 1, :]) / 4

    # Remainder path: 5 slices in groups of 2 → output of size 3 (2 full + 1).
    B = reshape(collect(1.0:30), 2, 3, 5)
    out_rem = group_z_project(B, 2; type = :sum)
    @test size(out_rem) == (2, 3, 3)
    @test out_rem[1, 1, 3] == B[1, 1, 5]

    # Max / min / prod.
    @test group_z_project(A, 4; type = :max)[1, 1, 1] == maximum(A[1, 1, :])
    @test group_z_project(A, 4; type = :min)[1, 1, 1] == minimum(A[1, 1, :])
    @test group_z_project(A, 2; type = :prod)[1, 1, 1] == A[1, 1, 1] * A[1, 1, 2]

    # Unsupported type
    @test_throws ArgumentError group_z_project(A, 2; type = :nope)
end

@testset "width_estimate — analytic FWHM of a Gaussian" begin
    σ = 5.0
    xs = -50:50
    v = exp.(-(xs .^ 2) ./ (2 * σ^2))
    fwhm_analytic = 2 * σ * sqrt(2 * log(2))
    @test isapprox(width_estimate(v, 0.5), fwhm_analytic; rtol = 0.02)

    # 1/e² width is wider than FWHM.
    @test width_estimate(v, exp(-2)) > width_estimate(v, 0.5)
end

@testset "width_estimate_3d — separable Gaussian" begin
    xs = -25:25
    σx, σy, σz = 3.0, 5.0, 7.0
    gx = exp.(-(xs .^ 2) ./ (2 * σx^2))
    gy = exp.(-(xs .^ 2) ./ (2 * σy^2))
    gz = exp.(-(xs .^ 2) ./ (2 * σz^2))
    M = [gx[i] * gy[j] * gz[k] for i in 1:51, j in 1:51, k in 1:51]
    wx, wy, wz = width_estimate_3d(M, 0.5)
    @test isapprox(wx, 2 * σx * sqrt(2 * log(2)); rtol = 0.05)
    @test isapprox(wy, 2 * σy * sqrt(2 * log(2)); rtol = 0.05)
    @test isapprox(wz, 2 * σz * sqrt(2 * log(2)); rtol = 0.05)
end

@testset "tpm_signal_scale — Xu–Webb formula" begin
    tpm = NAOMi.finalize!(NAOMi.TPMParams())
    psf = NAOMi.PSFParams()
    Ft = tpm_signal_scale(tpm, psf)
    @test Ft > 0
    @test isfinite(Ft)

    # Power-squared scaling.
    tpm2 = NAOMi.finalize!(NAOMi.TPMParams(pavg = tpm.pavg * 2))
    Ft2 = tpm_signal_scale(tpm2, psf)
    @test isapprox(Ft2 / Ft, 4.0; rtol = 1e-10)

    # No psf_params → reads from tpm.{nidx, nac, lambda}.
    @test tpm_signal_scale(tpm) > 0
end

@testset "collection_mask — hemoabs absorption" begin
    vp  = NAOMi.finalize!(NAOMi.VolumeParams(vol_sz = [30.0, 30.0, 10.0],
                                              vres = 2.0, vol_depth = 20.0))
    psf = NAOMi.PSFParams()
    volpx_xy = round(Int, vp.vol_sz[1] * vp.vres)
    volpx_z  = round(Int, (vp.vol_sz[3] + vp.vol_depth) * vp.vres)

    # Empty vasculature → unit transmission everywhere.
    ves_empty = zeros(Float64, volpx_xy, volpx_xy, volpx_z)
    cm0 = collection_mask(vp, psf, ves_empty)
    @test size(cm0) == (volpx_xy, volpx_xy)
    @test all(cm0 .≈ 1.0)

    # Uniform vasculature → uniform attenuation in the interior. The expected
    # central-pixel transmission is `10^(-N_eff / vres * hemoabs)`, where
    # `N_eff` is the column-integrated vessel path that falls inside the
    # collection cone. Verify it is less than 1 everywhere and that
    # deeper-depth volumes attenuate more.
    ves_full = ones(Float64, volpx_xy, volpx_xy, volpx_z)
    cm1 = collection_mask(vp, psf, ves_full)
    @test all(cm1 .< 1.0)
    @test all(cm1 .> 0.0)

    vp2 = NAOMi.finalize!(NAOMi.VolumeParams(vol_sz = [30.0, 30.0, 10.0],
                                              vres = 2.0, vol_depth = 40.0))
    volpx_z2 = round(Int, (vp2.vol_sz[3] + vp2.vol_depth) * vp2.vres)
    ves_full2 = ones(Float64, volpx_xy, volpx_xy, volpx_z2)
    cm2 = collection_mask(vp2, psf, ves_full2)
    @test minimum(cm2) < minimum(cm1)
end

@testset "collection_mask — hemoabs scaling exactness" begin
    # With a uniform-1 vasculature, the post-convolution column density at an
    # interior pixel equals the unfiltered z-sum (normalised disk kernel * 1
    # = 1 in the interior). So colmask(interior) = sum(coldist > 0
    # contributions) * proppx, and the result must equal
    # 10^(−col_density · hemoabs / vres).
    vp  = NAOMi.finalize!(NAOMi.VolumeParams(vol_sz = [30.0, 30.0, 10.0],
                                              vres = 2.0, vol_depth = 20.0))
    psf = NAOMi.PSFParams()
    volpx_xy = round(Int, vp.vol_sz[1] * vp.vres)
    volpx_z  = round(Int, (vp.vol_sz[3] + vp.vol_depth) * vp.vres)
    ves_full = ones(Float64, volpx_xy, volpx_xy, volpx_z)
    cm = collection_mask(vp, psf, ves_full)
    # Predict: proppx=20, blocks: 3, coldist (in pixels) values at iz=1,2,3
    proppx = round(Int, psf.prop_sz * vp.vres)
    cone_slope = vp.vres * tan(asin(psf.objNA / psf.n))
    coldist = [cone_slope * (vp.vol_depth + vp.vol_sz[3] / 2 -
                              proppx / vp.vres * (iz - 0.5))
               for iz in 1:cld(volpx_z, proppx)]
    expected_density = sum(proppx for c in coldist if c > 0)
    expected_transmission = 10.0^(-expected_density / vp.vres * psf.hemoabs)
    @test isapprox(cm[volpx_xy ÷ 2 + 1, volpx_xy ÷ 2 + 1],
                   expected_transmission; rtol = 1e-10)
end
