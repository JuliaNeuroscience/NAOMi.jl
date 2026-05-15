using NAOMi
using Test
using Random
using Statistics

@testset "samp_small_world_mat — basic shape" begin
    rng = MersenneTwister(42)
    A = samp_small_world_mat(rng, 30, 10, 0.0; self_ex = 0, rand_opt = 0)
    @test size(A) == (30, 30)
    # Toeplitz lattice from a length-`K_conn/2` block of ones — interior
    # rows have `2·(K_conn/2) − 1 = K_conn − 1` ones (the diagonal is
    # double-counted in the Toeplitz construction).
    interior = sum(A; dims = 2)[5:26]
    @test all(interior .== 9)
    # With self_ex added, the diagonal picks up the offset (lattice already
    # carries a 1.0 on the diagonal).
    A2 = samp_small_world_mat(rng, 30, 10, 0.0; self_ex = 4, rand_opt = 0)
    @test all([A2[i, i] for i in 1:30] .== 5.0)
end

@testset "samp_small_world_mat — N_node tuple appends background" begin
    rng = MersenneTwister(1)
    A = samp_small_world_mat(rng, (20, 5), 10, 0.0; self_ex = 0, rand_opt = 0)
    @test size(A) == (25, 25)
    # bg block (rows 21:25) is full of ones in the soma columns
    @test all(A[21:25, 1:20] .== 1.0)
    # bg-to-bg block is identity
    @test A[21:25, 21:25] == [i == j ? 1.0 : 0.0 for i in 1:5, j in 1:5]
end

@testset "samp_small_world_mat — argument validation" begin
    @test_throws ArgumentError samp_small_world_mat(10, 5, 0.0)
    @test_throws ArgumentError samp_small_world_mat(4, 10, 0.0)
end

@testset "samp_small_world_mat — spatial seeding picks nearest neighbours" begin
    rng = MersenneTwister(7)
    n_locs = reshape(collect(1.0:30.0), 30, 1)         # 1-D ordered locations
    A = samp_small_world_mat(rng, 30, 4, 0.0; self_ex = 0, rand_opt = 0,
                              n_locs = n_locs)
    # Each row should connect to itself + 3 nearest neighbours (i.e. ±1, ±2)
    # depending on edge effects; row 15 (well-interior) should connect to
    # {13, 14, 15, 16}.
    @test sort(findall(==(1.0), A[15, :])) ⊆ [13, 14, 15, 16, 17]
    @test sum(A[15, :]) == 4
end

@testset "expression_variation — scalar and 2-vector min_mod" begin
    rng = MersenneTwister(0)
    # Scalar form: uniform on [min_mod, 1] (when not silenced)
    x = expression_variation(rng, 1000, 0.0, 0.5)
    @test length(x) == 1000
    @test all(x .≥ 0.5)
    @test all(x .≤ 1.0)
    # p_off silences cells; with p_off = 0.9, expect ~10% non-zero
    rng = MersenneTwister(1)
    x = expression_variation(rng, 5000, 0.9, 0.5)
    @test count(==(0.0), x) > 4000
    # 2-vector form: Gamma-distributed (always positive when not silenced)
    rng = MersenneTwister(2)
    x = expression_variation(rng, 2000, 0.0, [0.4, 2.53])
    @test all(x .≥ 0)
    @test isapprox(mean(x), 0.4 * 2.53; rtol = 0.1)        # Gamma mean = α·θ
    # Argument validation
    @test_throws ArgumentError expression_variation(10, -0.1, 0.5)
    @test_throws ArgumentError expression_variation(10, 0.0, 1.5)
end

@testset "gen_correlated_spike_trains — discrete Hawkes" begin
    rng = MersenneTwister(0)
    so = SpikeOptions(K = 30, nt = 5000, dt = 1 / 100, rate = 1e-1,
                       burst_mean = 10.0, N_bg = 0)
    res = gen_correlated_spike_trains(rng, so)
    @test size(res.soma) == (30, 5000)
    @test size(res.bg)   == (0, 5000)
    @test size(res.A)    == (30, 30)
    @test sum(res.soma) > 0
    # Mean firing rate should be > 0; tight bound is hard for Hawkes.
    # Sanity: rates roughly in the range [0.5e-3, 5e-2] spikes/sample.
    emp = sum(res.soma) / length(res.soma)
    @test 1e-4 < emp < 1e-1
end

@testset "gen_correlated_spike_trains — N_bg > 0 path" begin
    rng = MersenneTwister(3)
    so = SpikeOptions(K = 20, nt = 1000, dt = 1 / 100, rate = 1e-2,
                       burst_mean = 5.0, N_bg = 5)
    res = gen_correlated_spike_trains(rng, so)
    @test size(res.soma) == (20, 1000)
    @test size(res.bg)   == (5, 1000)
end

@testset "gen_correlated_spike_trains — spatial correlation" begin
    # 1-D placement makes 'near' / 'far' unambiguous. With small-world
    # spatial seeding (beta = 0.3), most connections of each neuron are
    # to its nearest spatial neighbours, so pairs within a few units
    # should co-fire more than distant pairs. Sample-level binary spikes
    # are too sparse for direct correlation; bin into 50-sample windows
    # (~0.5 s) to capture burst-coincidence.
    rng = MersenneTwister(11)
    K = 30
    n_locs = reshape(collect(1.0:K), K, 1)
    so = SpikeOptions(K = K, nt = 20000, dt = 1 / 100, rate = 2e-1,
                       burst_mean = 20.0, N_bg = 0)
    res = gen_correlated_spike_trains(rng, so; n_locs = n_locs)
    bin = 50
    nb = size(res.soma, 2) ÷ bin
    Sb = zeros(Float64, K, nb)
    for j in 1:nb
        Sb[:, j] = sum(res.soma[:, (j - 1) * bin + 1:j * bin]; dims = 2)
    end
    Sd = Sb .- mean(Sb; dims = 2)
    norms = sqrt.(sum(Sd .^ 2; dims = 2))
    near_corr = Float64[]
    far_corr  = Float64[]
    for i in 1:K, j in (i + 1):K
        n = norms[i] * norms[j]
        n == 0 && continue
        c = (Sd[i, :]' * Sd[j, :]) / n
        if abs(i - j) ≤ 2
            push!(near_corr, c)
        elseif abs(i - j) ≥ 15
            push!(far_corr, c)
        end
    end
    @test !isempty(near_corr) && !isempty(far_corr)
    @test mean(near_corr) > mean(far_corr)
end

@testset "generate_time_traces — Ca_DE end-to-end" begin
    rng = MersenneTwister(0)
    so = SpikeOptions(K = 8, nt = 300, dt = 1 / 30, rate = 5e-2,
                       burst_mean = 10.0, smod_flag = :independent,
                       dyn_type = :Ca_DE, dendflag = true, axonflag = false,
                       N_bg = 0, p_off = 0.0, min_mod = [0.4, 2.53])
    res = generate_time_traces(rng, so)
    @test size(res.soma) == (8, 300)
    @test res.dend !== nothing && size(res.dend) == (8, 300)
    @test res.bg === nothing
    @test length(res.mod_vals) == 8
    @test all(isfinite, res.soma)
    @test all(isfinite, res.dend)
end

@testset "generate_time_traces — AR1 / AR2 branches" begin
    rng = MersenneTwister(0)
    for dt in (:AR1, :AR2)
        so = SpikeOptions(K = 6, nt = 200, dt = 1 / 30, rate = 5e-2,
                           burst_mean = 0.0, smod_flag = :independent,
                           dyn_type = dt, dendflag = false, axonflag = false,
                           N_bg = 0, p_off = 0.0, min_mod = [1.0])
        res = generate_time_traces(rng, so)
        @test size(res.soma) == (6, 200)
        @test all(isfinite, res.soma)
    end
end

@testset "generate_time_traces — Hawkes + N_bg > 0" begin
    rng = MersenneTwister(0)
    so = SpikeOptions(K = 12, nt = 200, dt = 1 / 30, rate = 5e-2,
                       burst_mean = 10.0, smod_flag = :hawkes,
                       dyn_type = :Ca_DE, dendflag = false, axonflag = false,
                       N_bg = 3, p_off = 0.0, min_mod = [1.0])
    res = generate_time_traces(rng, so)
    @test size(res.soma) == (12, 200)
    @test res.bg !== nothing && size(res.bg, 1) == 3
end

@testset "generate_time_traces — mod_vals override" begin
    rng = MersenneTwister(0)
    K = 5
    fixed = [0.0, 1.0, 0.5, 1.0, 0.2]
    so = SpikeOptions(K = K, nt = 100, dt = 1 / 30, rate = 5e-2,
                       burst_mean = 5.0, smod_flag = :independent,
                       dyn_type = :Ca_DE, dendflag = false, axonflag = false,
                       N_bg = 0, p_off = 0.0, min_mod = [1.0])
    res = generate_time_traces(rng, so; mod_vals = fixed)
    # Row 1 should be identically zero (mod_vals[1] = 0)
    @test all(res.soma[1, :] .== 0)
    @test res.mod_vals == fixed
end
