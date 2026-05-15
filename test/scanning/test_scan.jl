using NAOMi
using Test
using Random
using Statistics

# Shared volume + PSF + activity for all tests in this file.
function _scan_setup(; seed=0, nt=3, motion=true)
    vol = VolumeParams(vol_sz=[30, 30, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2, N_bg=10)
    finalize!(vol)
    np = NeuronParams()
    vp = VasculatureParams(sourceFreq=200.0, vesFreq=[60.0, 80.0, 30.0])
    dp = DendriteParams()
    ap = AxonParams(maxlength=40.0, padsize=4, numbranches=1)
    bp = BackgroundParams(maxlength=40.0)
    rng = MersenneTwister(seed)
    nv = simulate_neural_volume(vol, np, vp, dp, ap, bp; rng=rng)
    psf = zeros(Float32, 15, 15, 11)
    for k in 1:11, j in 1:15, i in 1:15
        psf[i, j, k] = exp(-((i - 8)^2 + (j - 8)^2) / 8 - (k - 6)^2 / 8)
    end
    psf ./= sum(psf)
    sp = ScanParams(scan_buff=4, motion=motion)
    neur_act = (soma=ones(Float32, vol.N_neur, nt),
                dend=ones(Float32, vol.N_neur, nt),
                bg=ones(Float32, ap.N_proc, nt))
    return (; nv, psf, sp, neur_act, vol, ap, rng)
end

@testset "scan_volume — output shape" begin
    s = _scan_setup(; nt=3, motion=false)
    mov = scan_volume(s.nv, s.psf, s.neur_act, s.sp; rng=MersenneTwister(0))
    @test ndims(mov) == 3
    @test eltype(mov) === Float32
    H_out = (size(s.nv.neur_vol, 1) - 2 * s.sp.scan_buff) ÷ s.sp.sfrac
    W_out = (size(s.nv.neur_vol, 2) - 2 * s.sp.scan_buff) ÷ s.sp.sfrac
    @test size(mov) == (H_out, W_out, 3)
end

@testset "scan_volume — motion off → frames identical" begin
    s = _scan_setup(; nt=2, motion=false)
    mov = scan_volume(s.nv, s.psf, s.neur_act, s.sp; rng=MersenneTwister(0))
    @test mov[:, :, 1] == mov[:, :, 2]
end

@testset "scan_volume — motion on → frames differ but bounded" begin
    s = _scan_setup(; nt=2, motion=true)
    mov = scan_volume(s.nv, s.psf, s.neur_act, s.sp; rng=MersenneTwister(0))
    @test mov[:, :, 1] != mov[:, :, 2]
    # Frame-to-frame difference should be a small fraction of frame mean.
    diff = maximum(abs.(mov[:, :, 1] .- mov[:, :, 2]))
    @test diff < maximum(mov[:, :, 1])
end

@testset "scan_volume — mean intensity stable across frames" begin
    s = _scan_setup(; nt=5, motion=false)
    nps = NoiseParams()
    tp = TPMParams(); finalize!(tp)
    so = SpikeOptions()
    mov = scan_volume(s.nv, s.psf, s.neur_act, s.sp;
                      noise_params=nps, tpm_params=tp, spike_opts=so,
                      rng=MersenneTwister(0))
    means = [mean(@view mov[:, :, k]) for k in 1:5]
    # Without motion or activity drift, per-frame means should be ~equal.
    @test (maximum(means) - minimum(means)) / mean(means) < 0.2
end

@testset "scan_volume — clean & motion outputs" begin
    s = _scan_setup(; nt=3, motion=true)
    mov, mov_clean, mot_hist =
        scan_volume(s.nv, s.psf, s.neur_act, s.sp;
                    rng=MersenneTwister(0),
                    return_clean=true, return_motion=true)
    @test size(mov_clean) == size(mov)
    @test size(mot_hist) == (3, 3)
    # z-row is the depth offset; should be near vol_depth/2.
    @test all(mot_hist[3, :] .>= 1)
end

@testset "img_sub_row_shift — zero offsets return centre crop" begin
    img = reshape(Float32.(1:100), 10, 10)
    out = img_sub_row_shift(img, 2, 0, zeros(Float32, 10))
    # buf_sz of 2 trims rows 1-2 and 9-10 and same in y → 6×6 output.
    @test size(out) == (10 - 4, 10 - 4)
end
