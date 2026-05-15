using NAOMi
using Test
using Random
using FFTW

@testset "psf_fft — round-trip via single_scan" begin
    # Delta-function PSF: FFT then single_scan should reproduce the input
    # (up to floating-point precision and the same-size crop).
    psf = zeros(Float32, 5, 5, 3)
    psf[3, 3, 2] = 1.0f0
    vol_sz = (10, 10, 5)
    fp = psf_fft(vol_sz, psf)
    @test size(fp, 3) == size(psf, 3)
    # FFT plane is zero-padded to vol_sz + psf_xy - 1.
    @test size(fp, 1) == vol_sz[1] + size(psf, 1) - 1
    @test size(fp, 2) == vol_sz[2] + size(psf, 2) - 1

    nv = zeros(Float32, vol_sz...)
    nv[5, 5, 2] = 1.0f0
    img = single_scan(nv, size(psf), fp; freq_opt=true)
    @test size(img) == (vol_sz[1], vol_sz[2])
    @test isapprox(sum(img), 1.0f0; atol=1e-5)
    @test argmax(img) == CartesianIndex(5, 5)
end

@testset "psf_fft + single_scan — z_sub presumming" begin
    psf = zeros(Float32, 3, 3, 4)
    psf[2, 2, :] .= 1.0f0
    fp = psf_fft((8, 8, 4), psf; z_sub=2)
    @test size(fp, 3) == 2
    nv = zeros(Float32, 8, 8, 4)
    nv[4, 4, 1] = 1.0f0; nv[4, 4, 3] = 2.0f0
    img = single_scan(nv, size(psf), fp; z_sub=2, freq_opt=true)
    # both source slices pair with PSF slices summing to 2 each → total 6.
    @test isapprox(sum(img), 6.0f0; atol=1e-4)
end

@testset "scan_volume_frame — end-to-end smoke" begin
    vol = VolumeParams(vol_sz=[30, 30, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2, N_bg=10)
    finalize!(vol)
    np = NeuronParams()
    vp = VasculatureParams(sourceFreq=200.0, vesFreq=[60.0, 80.0, 30.0])
    dp = DendriteParams()
    ap = AxonParams(maxlength=40.0, padsize=4, numbranches=1)
    bp = BackgroundParams(maxlength=40.0)
    rng = MersenneTwister(0)
    nv_vol = simulate_neural_volume(vol, np, vp, dp, ap, bp; rng=rng)

    psf = zeros(Float32, 15, 15, 11)
    σxy, σz = 2.0, 2.0
    for k in 1:11, j in 1:15, i in 1:15
        psf[i, j, k] = exp(-((i-8)^2 + (j-8)^2) / (2σxy^2) - (k-6)^2 / (2σz^2))
    end
    psf ./= sum(psf)

    sp = ScanParams(motion=false)
    scan_vol = setup_scan_volume_frame(nv_vol, psf, sp)
    neur_act = (soma=ones(Float32, vol.N_neur),
                dend=ones(Float32, vol.N_neur),
                bg=ones(Float32, length(scan_vol.axon_loc)))
    img, _ = scan_volume_frame(scan_vol, neur_act, sp)
    @test eltype(img) === Float32
    @test sum(img) > 0
    @test size(img, 1) == size(nv_vol.neur_vol, 1) ÷ sp.sfrac
    @test size(img, 2) == size(nv_vol.neur_vol, 2) ÷ sp.sfrac

    # Reproducibility: two calls on the same scan_vol give the same image.
    img2, _ = scan_volume_frame(scan_vol, neur_act, sp)
    @test img == img2
end

@testset "scan_volume_frame — TPM signal scaling (∝ pavg²)" begin
    vol = VolumeParams(vol_sz=[30, 30, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2, N_bg=5)
    finalize!(vol)
    np = NeuronParams()
    vp = VasculatureParams(sourceFreq=200.0, vesFreq=[60.0, 80.0, 30.0])
    dp = DendriteParams()
    ap = AxonParams(maxlength=40.0, padsize=4, numbranches=1)
    bp = BackgroundParams(maxlength=40.0)
    rng = MersenneTwister(0)
    nv_vol = simulate_neural_volume(vol, np, vp, dp, ap, bp; rng=rng)

    psf = zeros(Float32, 15, 15, 11)
    σxy = 2.0
    for k in 1:11, j in 1:15, i in 1:15
        psf[i, j, k] = exp(-((i-8)^2 + (j-8)^2) / (2σxy^2) - (k-6)^2 / 8)
    end
    psf ./= sum(psf)
    sp = ScanParams(motion=false)
    scan_vol = setup_scan_volume_frame(nv_vol, psf, sp)
    neur_act = (soma=ones(Float32, vol.N_neur),
                dend=ones(Float32, vol.N_neur),
                bg=ones(Float32, length(scan_vol.axon_loc)))

    tp1 = TPMParams(pavg=20.0); finalize!(tp1)
    tp2 = TPMParams(pavg=40.0); finalize!(tp2)
    img1, _ = scan_volume_frame(scan_vol, neur_act, sp; tpm_params=tp1)
    img2, _ = scan_volume_frame(scan_vol, neur_act, sp; tpm_params=tp2)
    @test sum(img1) > 0
    @test sum(img2) > 0
    @test isapprox(sum(img2) / sum(img1), (40 / 20)^2; atol=0.1)
end
