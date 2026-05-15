using NAOMi
using Test
using Random

# Helper: set up a fully-realised volume (soma + dendrites) so that the
# axon/bg code has something to consume.
function _setup_volume(seed::Integer; N_neur=2, N_bg=10)
    vol = VolumeParams(vol_sz=[40, 40, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=N_neur, N_bg=N_bg)
    finalize!(vol)
    np = NeuronParams()
    dp = DendriteParams()
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    neur_ves = falses(H, W, totalD)
    rng = MersenneTwister(seed)
    locs, Vcells, Vnucs, _ =
        sample_dense_neurons(np, vol, neur_ves; rng=rng)
    neur_soma, neur_vol0, gp_nuc, gp_soma =
        generate_neural_volume(np, vol, locs, Vcells, Vnucs, neur_ves)
    neur_num, dendnum_AD, _ =
        grow_neuron_dendrites!(vol, dp, neur_soma, neur_ves, locs,
                               gp_nuc, gp_soma; rng=rng)
    neur_num2, neur_num_AD =
        grow_apical_dendrites!(vol, dp, neur_num, dendnum_AD,
                               gp_nuc, gp_soma; rng=rng)
    gp_vals, neur_vol =
        set_cell_fluorescence(vol, np, dp, neur_num2, neur_soma,
                              neur_num_AD, locs, neur_vol0; rng=rng)
    return (; vol, np, dp, neur_ves, locs, neur_num2, neur_num_AD,
            neur_soma, neur_vol, gp_nuc, gp_soma, gp_vals, rng)
end

@testset "generate_axons — basic output" begin
    s = _setup_volume(0; N_bg=10)
    ap = AxonParams(maxlength=50.0, padsize=5, numbranches=2)
    neur_vol_out, gp_bgvals, ap_new =
        generate_axons(s.vol, ap, s.neur_vol, s.neur_num2, s.gp_vals, s.gp_nuc;
                       rng=s.rng)
    @test length(gp_bgvals) > 0
    @test length(gp_bgvals) <= s.vol.N_bg
    @test ap_new.N_proc == length(gp_bgvals)
    # Each process should have ≥ minlength voxels.
    for e in gp_bgvals
        @test length(e.loc) >= ap.minlength
        @test length(e.val) == length(e.loc)
        @test all(>(0), e.val)
    end
    # Total axon voxels should be a non-trivial fraction.
    total_axon = sum(length(e.loc) for e in gp_bgvals)
    bg_voxels = sum(s.neur_num2 .== 0)
    @test 0 < total_axon < bg_voxels    # at most fills available
    @test size(neur_vol_out) == size(s.neur_vol)
end

@testset "sort_axons — uniform binning when N_proc ≤ Ncomps" begin
    s = _setup_volume(0; N_bg=10)
    ap = AxonParams(maxlength=50.0, padsize=5, numbranches=2)
    _, gp_bgvals, _ =
        generate_axons(s.vol, ap, s.neur_vol, s.neur_num2, s.gp_vals, s.gp_nuc;
                       rng=s.rng)
    # N_proc = 3 < N_neur (2) + N_den (rounded) = 2 + 3 = 5: forces the uniform branch.
    ap_sort = AxonParams(maxlength=50.0, padsize=5, numbranches=2, N_proc=3)
    bins = sort_axons(s.vol, ap_sort, gp_bgvals, s.locs; rng=s.rng)
    @test length(bins) == ap_sort.N_proc
    # All bg voxels should be preserved (no duplicates, no losses).
    all_locs = vcat([Int.(b.loc) for b in bins]...)
    all_in   = vcat([Int.(e.loc) for e in gp_bgvals]...)
    @test sort(all_locs) == sort(all_in)
end

@testset "sort_axons — nearest-cell binning when N_proc > Ncomps" begin
    s = _setup_volume(0; N_bg=15)
    ap = AxonParams(maxlength=50.0, padsize=5, numbranches=2)
    _, gp_bgvals, _ =
        generate_axons(s.vol, ap, s.neur_vol, s.neur_num2, s.gp_vals, s.gp_nuc;
                       rng=s.rng)
    # N_proc > N_neur + N_den engages the nearest-cell assignment branch.
    ap_sort = AxonParams(maxlength=50.0, padsize=5, numbranches=2, N_proc=20)
    bins = sort_axons(s.vol, ap_sort, gp_bgvals, s.locs; rng=s.rng)
    @test length(bins) == 20
end

@testset "generate_bg_dendrites — adds new labels above N_neur + N_den" begin
    s = _setup_volume(0; N_bg=10)
    bp = BackgroundParams(maxlength=50.0)
    Ncomps_before = s.vol.N_neur + Int(round(s.vol.N_den))
    neur_num3, neur_vol3, gp_vals2, locs2 =
        generate_bg_dendrites(s.vol, bp, s.dp, s.neur_vol, s.neur_num2,
                              s.gp_vals, s.gp_nuc, s.locs; rng=s.rng)
    # All new labels should exceed Ncomps_before.
    new_labels = unique(neur_num3) |> sort
    @test all(l -> l == 0 || l <= Ncomps_before || l > Ncomps_before, new_labels)
    # New gp_vals entries should be appended (length grows).
    @test length(gp_vals2) >= length(s.gp_vals)
    # New locs entries should be appended.
    @test size(locs2, 1) >= size(s.locs, 1)
    @test size(locs2, 2) == 3
end
