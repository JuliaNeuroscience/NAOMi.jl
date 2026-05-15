using NAOMi
using Test
using Random

@testset "simulate_neural_volume — end-to-end smoke" begin
    vol = VolumeParams(vol_sz=[30, 30, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2, N_bg=10)
    finalize!(vol)
    np = NeuronParams()
    vp = VasculatureParams(sourceFreq=200.0, vesFreq=[60.0, 80.0, 30.0])
    dp = DendriteParams()
    ap = AxonParams(maxlength=40.0, padsize=4, numbranches=1)
    bp = BackgroundParams(maxlength=40.0)
    rng = MersenneTwister(0)

    nv = simulate_neural_volume(vol, np, vp, dp, ap, bp; rng=rng)

    # Component-array dimensions are all consistent.
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    D = Int(vol.vol_sz[3] * vol.vres)
    @test size(nv.neur_vol)   == (H, W, D)
    @test size(nv.neur_num)   == (H, W, D)
    @test size(nv.neur_soma)  == (H, W, D)
    @test size(nv.neur_num_AD) == (H, W, D)
    @test size(nv.neur_ves)   == (H, W, D)

    # Output types.
    @test eltype(nv.neur_vol)    === Float32
    @test eltype(nv.neur_num)    === UInt16
    @test eltype(nv.neur_soma)   === UInt16
    @test eltype(nv.neur_num_AD) === UInt16

    # Reasonable population.
    @test sum(nv.neur_soma .> 0) > 0
    @test sum(nv.neur_num .> 0)  > sum(nv.neur_soma .> 0)
    @test all(0 .<= nv.neur_vol)
    # Soma and dendrite-only voxel sets do not overlap inside neur_soma.
    @test count(nv.neur_soma .> 0) == count(nv.neur_soma .> 0 .& (nv.neur_num .> 0))

    # gp_vals length equals N_neur + N_den (+ any bg-dendrite extras).
    @test length(nv.gp_vals) >= vol.N_neur + Int(round(vol.N_den))
    # gp_nuc carries N_neur entries.
    @test length(nv.gp_nuc) == vol.N_neur

    # bg_proc count.
    @test length(nv.bg_proc) == ap.N_proc

    # Soma centres (first N_neur rows) are in the brain-volume slab; bg
    # dendrite roots (later rows) are placed *outside* the volume.
    soma_locs = nv.locs[1:vol.N_neur, :]
    @test all(0 .<= soma_locs[:, 1] .<= vol.vol_sz[1] + 0.5)
    @test all(0 .<= soma_locs[:, 2] .<= vol.vol_sz[2] + 0.5)
    @test all(0 .<= soma_locs[:, 3] .<= vol.vol_sz[3] + 0.5)
end

@testset "simulate_neural_volume — flags disable subsystems" begin
    vol = VolumeParams(vol_sz=[30, 30, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2, N_bg=5)
    finalize!(vol)
    np = NeuronParams()
    vp = VasculatureParams(flag=0)
    dp = DendriteParams()
    ap = AxonParams(flag=0)
    bp = BackgroundParams(flag=0)
    rng = MersenneTwister(0)

    nv = simulate_neural_volume(vol, np, vp, dp, ap, bp; rng=rng)

    # Vasculature off → no vessel voxels.
    @test sum(nv.neur_ves) == 0
    # Axons off → empty gp_bgvals / bg_proc.
    @test isempty(nv.gp_bgvals)
    @test isempty(nv.bg_proc)
    # Somas + dendrites still present.
    @test sum(nv.neur_soma .> 0) > 0
end
