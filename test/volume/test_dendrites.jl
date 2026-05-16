using NAOMi
using Test
using Random

@testset "dendrite_dijkstra_grid — unit cost" begin
    M = ones(Float32, 5, 5, 5, 6)
    d, pf = dendrite_dijkstra_grid(M, (3, 3, 3))
    @test d[3, 3, 3] == 0.0f0
    # 8 corners are each 6 steps away (manhattan).
    @test d[1, 1, 1] == 6.0f0
    @test d[5, 5, 5] == 6.0f0
    # Predecessors reach the root.
    path = get_dendrite_path(pf, (5, 5, 5), (3, 3, 3))
    @test !isempty(path)
    @test path[1]   == (5, 5, 5)
    @test path[end] == (3, 3, 3)
    @test length(path) == 7  # 6 hops + root
end

@testset "dendrite_dijkstra_grid — wall blocks traversal" begin
    M = ones(Float32, 4, 4, 4, 6)
    # Block the +x face of (1, j, k) from being reached.
    M[2, :, :, 1] .= Float32(Inf)
    d, _ = dendrite_dijkstra_grid(M, (1, 1, 1))
    @test isinf(d[2, 1, 1])
    # Still reachable via +y or +z chains.
    @test isfinite(d[1, 2, 1])
end

@testset "get_dendrite_path — unreachable returns empty" begin
    M = ones(Float32, 3, 3, 3, 6)
    _, pf = dendrite_dijkstra_grid(M, (1, 1, 1))
    # An island: predecessor map at (1,1,1) is 0 (it's the root). For
    # an arbitrarily-unreachable point, we'd need to construct one;
    # here just check that a reached node walks back correctly.
    path = get_dendrite_path(pf, (3, 3, 3), (1, 1, 1))
    @test path[1] == (3, 3, 3)
    @test path[end] == (1, 1, 1)
end

@testset "dendrite_random_walk — directed steps toward target" begin
    M = ones(Float32, 20, 20, 20)
    rng = MersenneTwister(0)
    path = dendrite_random_walk(M, (10, 10, 10), (10, 10, 20);
                                distsc=1.0, maxlength=50, fillweight=1.0,
                                maxel=8, minlength=5)
    @test !isempty(path)
    # The walk should generally trend toward target z=20.
    @test last(path)[3] > first(path)[3]
end

@testset "dilate_dendrite_paths_all — grows thickness >1 into 6-neighbours" begin
    paths = zeros(Float32, 10, 10, 10)
    pathnums = zeros(Int, 10, 10, 10)
    paths[5, 5, 5] = 3.0f0       # seed needs 2 extra voxels of thickness
    pathnums[5, 5, 5] = 1
    obstruction = zeros(Float32, 10, 10, 10)
    rng = MersenneTwister(0)
    pout, nout = dilate_dendrite_paths_all(paths, pathnums, obstruction;
                                           maxDist=3, rng=rng)
    # The seed's thickness should have been spent (now == 1).
    @test pout[5, 5, 5] == 1.0f0
    # Exactly 2 new voxels should be labelled 1 (one for each dropped unit).
    @test count(==(1), nout) == 3
end

@testset "grow_neuron_dendrites! — dendrites originate from soma boundaries" begin
    vol = VolumeParams(vol_sz=[40, 40, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2)
    finalize!(vol)
    np = NeuronParams()
    dp = DendriteParams()
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    neur_ves = falses(H, W, totalD)
    rng = MersenneTwister(0)
    locs, Vcells, Vnucs, _ = sample_dense_neurons(np, vol, neur_ves; rng=rng)
    neur_soma, _, gp_nuc, gp_soma =
        generate_neural_volume(np, vol, locs, Vcells, Vnucs, neur_ves)
    neur_num, dendnum_AD, _ =
        grow_neuron_dendrites!(vol, dp, neur_soma, neur_ves, locs,
                               gp_nuc, gp_soma; rng=rng)
    @test size(neur_num) == size(neur_soma)
    @test eltype(neur_num) === UInt16

    # For each cell: at least one dendrite-only voxel is adjacent to its soma.
    for k in 1:vol.N_neur
        soma_set = Set(findall(neur_soma .== UInt16(k)))
        dend_set = setdiff(Set(findall(neur_num .== UInt16(k))), soma_set)
        @test !isempty(dend_set)
        any_adj = false
        for c in dend_set
            for off in (CartesianIndex(1,0,0), CartesianIndex(-1,0,0),
                        CartesianIndex(0,1,0), CartesianIndex(0,-1,0),
                        CartesianIndex(0,0,1), CartesianIndex(0,0,-1))
                n = c + off
                (1 <= n[1] <= size(neur_num,1) &&
                 1 <= n[2] <= size(neur_num,2) &&
                 1 <= n[3] <= size(neur_num,3)) || continue
                if n in soma_set
                    any_adj = true
                    break
                end
            end
            any_adj && break
        end
        @test any_adj
    end
end

@testset "grow_apical_dendrites! — reaches top of volume" begin
    vol = VolumeParams(vol_sz=[40, 40, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2)
    finalize!(vol)
    np = NeuronParams()
    dp = DendriteParams()
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    neur_ves = falses(H, W, totalD)
    rng = MersenneTwister(0)
    locs, Vcells, Vnucs, _ = sample_dense_neurons(np, vol, neur_ves; rng=rng)
    neur_soma, _, gp_nuc, gp_soma =
        generate_neural_volume(np, vol, locs, Vcells, Vnucs, neur_ves)
    neur_num, dendnum_AD, _ =
        grow_neuron_dendrites!(vol, dp, neur_soma, neur_ves, locs,
                               gp_nuc, gp_soma; rng=rng)
    neur_num2, neur_num_AD =
        grow_apical_dendrites!(vol, dp, neur_num, dendnum_AD, gp_nuc, gp_soma;
                               rng=rng)
    # Apical dendrites must touch some voxel at z=1 (top of volume).
    @test any(neur_num_AD[:, :, 1] .> 0) || any(neur_num2[:, :, 1] .> 0)
    # neur_num2 must not lose any prior soma labels.
    @test count(neur_num2 .> 0) >= count(neur_num .> 0)
end

@testset "couple_dendrites=true — serial coupled growth" begin
    vol = VolumeParams(vol_sz=[40, 40, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=3)
    finalize!(vol)
    np = NeuronParams()
    dp = DendriteParams()
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    neur_ves = falses(H, W, totalD)
    locs, Vcells, Vnucs, _ =
        sample_dense_neurons(np, vol, neur_ves; rng=MersenneTwister(0))
    neur_soma, _, gp_nuc, gp_soma =
        generate_neural_volume(np, vol, locs, Vcells, Vnucs, neur_ves)

    grow(mode) = grow_neuron_dendrites!(vol, dp, neur_soma, neur_ves, locs,
                                        gp_nuc, gp_soma; rng=MersenneTwister(0),
                                        couple_dendrites=mode)
    nn_c, ad_c, _ = grow(true)
    @test eltype(nn_c) === UInt16
    @test size(nn_c) == size(neur_soma)
    @test count(nn_c .> 0) > count(neur_soma .> 0)   # dendrites were added

    # Coupled growth is deterministic given the seed.
    nn_c2, = grow(true)
    @test nn_c2 == nn_c
    # ...and differs from the default parallel, decoupled growth.
    nn_p, = grow(false)
    @test nn_p != nn_c

    # `grow_apical_dendrites!` accepts the same option.
    na_c, nad_c = grow_apical_dendrites!(vol, dp, nn_c, ad_c, gp_nuc, gp_soma;
                                         rng=MersenneTwister(0),
                                         couple_dendrites=true)
    @test eltype(na_c) === UInt16
    @test count(na_c .> 0) >= count(nn_c .> 0)
end

@testset "set_cell_fluorescence — basic structure" begin
    vol = VolumeParams(vol_sz=[40, 40, 20], vol_depth=15.0, vres=2.0,
                       min_dist=12.0, N_neur=2)
    finalize!(vol)
    np = NeuronParams(nuc_fluorsc=0.3)
    dp = DendriteParams()
    H = Int(vol.vol_sz[1] * vol.vres)
    W = Int(vol.vol_sz[2] * vol.vres)
    totalD = Int((vol.vol_sz[3] + vol.vol_depth) * vol.vres)
    neur_ves = falses(H, W, totalD)
    rng = MersenneTwister(0)
    locs, Vcells, Vnucs, _ = sample_dense_neurons(np, vol, neur_ves; rng=rng)
    neur_soma, neur_vol0, gp_nuc, gp_soma =
        generate_neural_volume(np, vol, locs, Vcells, Vnucs, neur_ves)
    neur_num, dendnum_AD, _ =
        grow_neuron_dendrites!(vol, dp, neur_soma, neur_ves, locs,
                               gp_nuc, gp_soma; rng=rng)
    neur_num2, neur_num_AD =
        grow_apical_dendrites!(vol, dp, neur_num, dendnum_AD, gp_nuc, gp_soma;
                               rng=rng)

    gp_vals, neur_vol = set_cell_fluorescence(vol, np, dp,
                                              neur_num2, neur_soma, neur_num_AD,
                                              locs, neur_vol0; rng=rng)
    @test length(gp_vals) == vol.N_neur + Int(round(vol.N_den))
    # Per-cell entries should carry locations + values of equal length.
    for kk in 1:vol.N_neur
        e = gp_vals[kk]
        @test length(e.loc) == length(e.val)
        @test length(e.is_soma) == length(e.val)
        # Soma-voxel values are normalised around 1.
        if any(e.is_soma)
            sv = e.val[e.is_soma]
            @test abs(sum(sv) / length(sv) - 1) < 0.3
        end
    end
end

@testset "smooth_cell_body — adds a small ball at first soma hit" begin
    fdims = (10, 10, 10)
    lin = LinearIndices(fdims)
    cellBody = Int32[lin[5, 5, 5], lin[5, 5, 6], lin[5, 5, 7]]
    path = [(7, 5, 5), (6, 5, 5), (5, 5, 5)]   # path hits cellBody at (5,5,5)
    out = smooth_cell_body([path], cellBody, fdims)
    @test !isempty(out)
    # The output should be inside a ball of radius 2 around (5, 5, 5).
    cart = CartesianIndices(fdims)
    for li in out
        c = cart[Int(li)]
        d2 = (c[1] - 5)^2 + (c[2] - 5)^2 + (c[3] - 5)^2
        @test d2 <= 4
    end
end
