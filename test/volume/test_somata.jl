using NAOMi
using Test
using Random
using LinearAlgebra

@testset "spiral_sample_sphere" begin
    V = spiral_sample_sphere(200)
    @test size(V) == (200, 3)
    radii = vec(sqrt.(sum(V .^ 2; dims=2)))
    @test all(isapprox.(radii, 1.0; atol=1e-10))
    # Approximate uniformity: 8 equal-area z-bands should each carry the
    # same count to within √N.
    bin = clamp.(fld.(V[:, 3] .+ 1.0, 2 / 8) .+ 1, 1, 8)
    counts = [count(==(b), bin) for b in 1:8]
    @test maximum(counts) - minimum(counts) <= ceil(sqrt(200))
end

@testset "teardrop_projection" begin
    V = spiral_sample_sphere(100)
    T = teardrop_projection(V, 1)
    @test size(T) == size(V)
    # z = -cos(θ) maps the south pole to +1 and north pole to -1.
    @test isapprox(T[end, 3],  1.0; atol=1e-6)
    @test isapprox(T[1,   3], -1.0; atol=1e-6)
    # NaN scrubbing: no NaNs anywhere.
    @test !any(isnan, T)
end

@testset "generate_neural_body" begin
    np = NeuronParams()
    rng = MersenneTwister(0)
    Vc, Vn, rot = generate_neural_body(np; rng=rng)
    @test size(Vc) == (np.n_samps, 3)
    @test size(Vn) == (np.n_samps, 3)
    @test length(rot) == 3
    @test all(abs.(rot) .<= np.max_ang + 1e-9)

    radii_c = vec(sqrt.(sum(Vc .^ 2; dims=2)))
    @test minimum(radii_c) > 0
    # Within the upstream-prescribed soma-radius bracket plus
    # nuc-offset slack.
    @test minimum(radii_c) >= np.exts[1] * np.avg_rad - 5.0
    @test maximum(radii_c) <= np.exts[2] * np.avg_rad + 5.0

    # Per-axis bounding-box ratio is bounded by eccens × exts (see soma
    # generation notes); a generous ratio cap catches degenerate shapes.
    bbox = maximum(Vc; dims=1) .- minimum(Vc; dims=1)
    ratio = maximum(bbox) / minimum(bbox)
    @test ratio < 5.0
end

@testset "sample_dense_neurons — count + min-distance" begin
    vol = VolumeParams(vol_sz=[60, 60, 30], vol_depth=20.0, vres=2.0,
                       min_dist=16.0, N_neur=10)
    finalize!(vol)
    np  = NeuronParams()
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    neur_ves = falses(H, W, totalD)
    rng = MersenneTwister(0)
    locs, Vcells, Vnucs, rots = sample_dense_neurons(np, vol, neur_ves; rng=rng)
    @test size(locs) == (vol.N_neur, 3)
    @test length(Vcells) == vol.N_neur
    @test length(Vnucs) == vol.N_neur
    # Minimum pairwise distance respected.
    mind = Inf
    for i in 1:size(locs, 1), j in (i + 1):size(locs, 1)
        d = sqrt(sum((locs[i, :] - locs[j, :]) .^ 2))
        mind = min(mind, d)
    end
    @test mind >= vol.min_dist
end

@testset "sample_dense_neurons — vasculature avoidance" begin
    vol = VolumeParams(vol_sz=[40, 40, 20], vol_depth=10.0, vres=2.0,
                       min_dist=12.0, N_neur=5)
    finalize!(vol)
    np = NeuronParams()
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    # Block off the bottom half of the brain volume with vasculature.
    neur_ves = falses(H, W, totalD)
    vol_depth_vox = Int(vol.vol_depth * vol.vres)
    neur_ves[:, :, (vol_depth_vox + Int(vol.vol_sz[3] * vol.vres / 2)):end] .= true
    rng = MersenneTwister(0)
    locs, _, _, _ = sample_dense_neurons(np, vol, neur_ves; rng=rng)
    # Every accepted soma should sit in the unblocked top half of the brain.
    @test all(locs[:, 3] .<= vol.vol_sz[3] / 2 + 1.0)
end

@testset "generate_neural_volume — rasterization" begin
    vol = VolumeParams(vol_sz=[60, 60, 30], vol_depth=20.0, vres=2.0,
                       min_dist=16.0, N_neur=5)
    finalize!(vol)
    np = NeuronParams(nuc_fluorsc=0.3)
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    neur_ves = falses(H, W, totalD)
    rng = MersenneTwister(0)
    locs, Vcells, Vnucs, _ = sample_dense_neurons(np, vol, neur_ves; rng=rng)
    neur_soma, neur_vol, gp_nuc, gp_soma =
        generate_neural_volume(np, vol, locs, Vcells, Vnucs, neur_ves)

    @test size(neur_soma) == (H, W, Int(vol.vol_sz[3] * vol.vres))
    @test eltype(neur_soma) === UInt16
    @test maximum(neur_soma) == UInt16(length(Vcells))
    @test sum(neur_soma .> 0) > 0
    # Each cell labelled in neur_soma should match the per-cell gp_soma list.
    for k in 1:length(gp_soma)
        labelled = count(==(UInt16(k)), neur_soma)
        @test labelled == length(gp_soma[k])
    end
    # Nucleus fluorescence ratio: gp_nuc[k][2] equals neur_params.nuc_fluorsc.
    for (idxs, val) in gp_nuc
        @test val == np.nuc_fluorsc
        if !isempty(idxs)
            @test all(neur_vol[i] == Float32(np.nuc_fluorsc) for i in idxs)
        end
    end
end

@testset "point_in_soma — star-shape radial test" begin
    # Synthetic mesh: 12 directions on unit sphere at radius 5, centered at origin
    V = spiral_sample_sphere(12) .* 5.0
    center = [0.0, 0.0, 0.0]
    @test point_in_soma(V, [0.0, 0.0, 0.0], center)
    @test point_in_soma(V, [2.0, 1.0, 0.0], center)
    @test !point_in_soma(V, [10.0, 0.0, 0.0], center)
    # off-center
    @test point_in_soma(V .+ [1.0 0.0 0.0], [3.0, 0.0, 0.0], [1.0, 0.0, 0.0])
end

@testset "isolate_visible_somas" begin
    # Gaussian PSF for a clean half-width.
    psf = zeros(Float32, 21, 21, 31)
    σxy, σz = 2.0, 4.0
    for k in 1:31, j in 1:21, i in 1:21
        psf[i, j, k] = exp(-((i - 11)^2 + (j - 11)^2) / (2σxy^2) - (k - 16)^2 / (2σz^2))
    end
    vol = VolumeParams(vol_sz=[60, 60, 30], vres=2.0, N_neur=10)
    finalize!(vol)
    np = NeuronParams(avg_rad=1.0)
    # Three somas: near midplane, near top, near bottom.
    locs = Float64[10 10 15; 20 20 5; 30 30 25]
    vis = isolate_visible_somas(locs, psf, vol, np)
    # The midplane one (z=15 == vol_sz[3]/2) is always within threshold.
    @test 1 in vis
end

@testset "masked_3d_gp" begin
    rng = MersenneTwister(0)
    gp = masked_3d_gp(8, 0.1, 1.0, 0.0; rng=rng)
    @test size(gp) == (8, 8, 8)
    @test eltype(gp) === Float32
    @test all(isfinite, gp)
    # Adding a positive mean shifts the field.
    gp2 = masked_3d_gp(8, 0.1, 1.0, 5.0; rng=MersenneTwister(0))
    @test all(isapprox.(gp2 .- gp, 5.0f0; atol=1e-4))
end
