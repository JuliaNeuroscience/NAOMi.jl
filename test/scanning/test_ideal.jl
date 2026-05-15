using NAOMi
using Test
using Random
using Statistics

# Build a small synthetic test volume, a Gaussian PSF, and a known
# soma-activity time course; verify that ideal profiles + times_from_profs
# recover the activity to high correlation.
function _ideal_setup(; seed=0, nt=40)
    vol = VolumeParams(vol_sz=[30, 30, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2, N_bg=2)
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

    sp = ScanParams(scan_buff=4, motion=false)
    K = vol.N_neur
    # Two soma traces: one slow sinusoid, one fast sinusoid, both with a
    # strictly positive minimum so that calculate_ideal_comps gets a
    # non-zero baseline activation for each cell.
    soma_act = zeros(Float32, K, nt)
    for t in 1:nt
        soma_act[1, t] = 1.0f0 + 0.5f0 * sin(2π * (t - 1) / nt)
        if K >= 2
            soma_act[2, t] = 1.0f0 + 0.5f0 * sin(4π * (t - 1) / nt + 1)
        end
    end
    dend_act = copy(soma_act)
    bg_act = fill(0.5f0, max(ap.N_proc, 1), nt)
    neur_act = (soma=soma_act, dend=dend_act, bg=bg_act)
    return (; nv, psf, sp, neur_act, vol, soma_act)
end

@testset "calculate_ideal_comps — shape and non-zero output" begin
    s = _ideal_setup(; nt=20)
    comps, baseim, ideal =
        calculate_ideal_comps(s.nv, s.psf, s.neur_act, s.sp;
                              rng=MersenneTwister(0))
    H_out = (size(s.nv.neur_vol, 1) - 2 * s.sp.scan_buff) ÷ s.sp.sfrac
    W_out = (size(s.nv.neur_vol, 2) - 2 * s.sp.scan_buff) ÷ s.sp.sfrac
    K = s.vol.N_neur
    @test size(comps) == (H_out, W_out, K)
    @test size(ideal) == (H_out, W_out, K)
    @test size(baseim) == (H_out, W_out)
    # FFT-based convolution leaves O(1e-8) negative roundoff in clean
    # output; assert against that floor rather than strict non-negativity.
    @test minimum(comps) > -1f-6
    @test sum(comps) > 0   # at least one cell produces a non-empty profile
end

@testset "comps2ideals — masks non-soma pixels to zero" begin
    H, W, K = 12, 12, 2
    comps = zeros(Float32, H, W, K)
    # Cell 1 has a clean 4×4 hot spot; cell 2 has a smaller one offset.
    comps[3:6, 3:6, 1] .= 1f0
    comps[8:10, 8:10, 2] .= 1f0
    baseim = fill(0.5f0, H, W)
    ideal = comps2ideals(comps, baseim)
    @test all(ideal .>= 0)
    # Each component slice retains only the hot spot.
    @test sum(ideal[3:6, 3:6, 1] .> 0) >= 5
    @test sum(ideal[:, :, 1] .> 0) == sum(ideal[3:6, 3:6, 1] .> 0)
    @test sum(ideal[8:10, 8:10, 2] .> 0) >= 5
    @test sum(ideal[:, :, 2] .> 0) == sum(ideal[8:10, 8:10, 2] .> 0)
end

@testset "times_from_profs — recovers ground-truth activity" begin
    H, W, T = 16, 16, 30
    # Two disjoint Gaussian-blob profiles.
    profs = zeros(Float32, H, W, 2)
    for j in 1:W, i in 1:H
        profs[i, j, 1] = exp(-((i - 5)^2 + (j - 5)^2) / 4)
        profs[i, j, 2] = exp(-((i - 12)^2 + (j - 12)^2) / 4)
    end
    truth = zeros(Float32, 2, T)
    for t in 1:T
        truth[1, t] = 1f0 + sin(2π * (t - 1) / T)
        truth[2, t] = 1f0 + cos(2π * (t - 1) / T)
    end
    mov = zeros(Float32, H, W, T)
    for t in 1:T, j in 1:W, i in 1:H
        mov[i, j, t] = profs[i, j, 1] * truth[1, t] +
                       profs[i, j, 2] * truth[2, t]
    end
    x_est, _ = times_from_profs(mov, profs; lambda=0)
    @test size(x_est) == (2, T)
    # Each recovered trace correlates strongly with the corresponding truth.
    for k in 1:2
        ρ = cor(vec(x_est[k, :]), vec(truth[k, :]))
        @test ρ > 0.95
    end
end

@testset "ideal pipeline — synthetic recovery from clean movie" begin
    s = _ideal_setup(; nt=30)
    # Build a clean (noise-free) movie with the same activity.
    mov_clean = scan_volume(s.nv, s.psf, s.neur_act, s.sp;
                            rng=MersenneTwister(1))
    # Spatial profiles via calculate_ideal_comps.
    comps, baseim, ideal =
        calculate_ideal_comps(s.nv, s.psf, s.neur_act, s.sp;
                              rng=MersenneTwister(0))
    # Recover traces from the clean movie against the comps profiles.
    x_est, _ = times_from_profs(mov_clean, comps; lambda=0)
    @test size(x_est, 1) == s.vol.N_neur
    @test size(x_est, 2) == size(s.soma_act, 2)
    # The cell with a non-empty profile should correlate with its truth.
    for k in 1:s.vol.N_neur
        col_sum = sum(comps[:, :, k])
        col_sum > 0 || continue
        ρ = cor(vec(x_est[k, :]), vec(s.soma_act[k, :]))
        @test ρ > 0.5
    end
end

@testset "times_from_profs — unconstrained LS path" begin
    H, W, T = 10, 10, 12
    profs = zeros(Float32, H, W, 1)
    profs[3:6, 3:6, 1] .= 1f0
    truth = reshape(Float32.(1:T), 1, T)
    mov = zeros(Float32, H, W, T)
    for t in 1:T
        mov[:, :, t] .= profs[:, :, 1] .* truth[1, t]
    end
    x_est, _ = times_from_profs(mov, profs; lambda=0, nnls=false)
    @test cor(vec(x_est[1, :]), vec(truth[1, :])) > 0.99
end
