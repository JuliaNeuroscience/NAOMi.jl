# Top-level neural-volume orchestrator.
#
# Ported from the upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - simulate_neural_volume.m → simulate_neural_volume
#
# Note: `resampVolume.m` is an empty stub upstream and is not ported.
# All other upstream Chunk-12 helpers (`gennode`, `delnode`,
# `nodesToConn`, `genconn`, `connToVol`, `branchGrowNodes`) were
# already ported in Chunk 8 (`src/volume/vasculature.jl`) to keep the
# orchestrator thin.

export simulate_neural_volume, NeuralVolume

"""
    NeuralVolume

Result of [`simulate_neural_volume`](@ref). All voxel arrays are sized
`(vol_sz[1]*vres, vol_sz[2]*vres, vol_sz[3]*vres)`. Public fields
mirror upstream `vol_out`:

| Field          | Description |
|----------------|-------------|
| `neur_vol`     | per-voxel base fluorescence (Float32) |
| `neur_num`     | per-voxel cell / dendrite / background label (UInt16) |
| `neur_soma`    | per-voxel soma label (UInt16, 0 outside somas) |
| `neur_num_AD`  | per-voxel apical-dendrite label (UInt16) |
| `gp_nuc`       | per-cell `(loc, val)` tuples for nucleus voxels |
| `gp_soma`      | per-cell soma voxel-index lists |
| `gp_vals`      | per-cell `(loc, val, is_soma)` fluorescence tuples (incl. bg) |
| `gp_bgvals`    | per-process axon `(loc, val)` tuples |
| `bg_proc`      | sorted background-process bins (`N_proc` entries) |
| `neur_ves`     | vessel mask (BitArray{3} sized as the *brain-only* slab) |
| `neur_ves_all` | full vasculature mask incl. above-volume slab (BitArray{3}) |
| `locs`         | `K×3` neuron centre positions (microns) |
| `Vcells`       | per-cell soma surface meshes |
| `Vnucs`        | per-cell nucleus surface meshes |
"""
struct NeuralVolume
    neur_vol::Array{Float32, 3}
    neur_num::Array{UInt16, 3}
    neur_soma::Array{UInt16, 3}
    neur_num_AD::Array{UInt16, 3}
    gp_nuc::Vector{Tuple{Vector{Int32}, Float64}}
    gp_soma::Vector{Vector{Int32}}
    gp_vals::Vector
    gp_bgvals::Vector
    bg_proc::Vector
    neur_ves::BitArray{3}
    neur_ves_all::BitArray{3}
    locs::Matrix{Float64}
    Vcells::Vector{Matrix{Float64}}
    Vnucs::Vector{Matrix{Float64}}
end

"""
    simulate_neural_volume(vol_params, neur_params, vasc_params, dend_params,
                           axon_params, bg_params; rng=Random.default_rng())
        -> NeuralVolume

End-to-end neural-volume generation pipeline. Calls (in order):

1. [`simulate_blood_vessels`](@ref) for the vasculature mask.
2. [`sample_dense_neurons`](@ref) for soma positions + meshes.
3. [`generate_neural_volume`](@ref) to rasterize somata and nuclei.
4. [`grow_neuron_dendrites!`](@ref) to add basal + local-apical dendrites.
5. [`grow_apical_dendrites!`](@ref) to add through-volume apical dendrites.
6. [`set_cell_fluorescence`](@ref) to populate per-cell fluorescence values.
7. [`generate_bg_dendrites`](@ref) (if `bg_params.flag != 0`) to add
   neuropil-background processes.
8. [`generate_axons`](@ref) + [`sort_axons`](@ref) (if
   `axon_params.flag != 0`) to add the per-process axon channel.

Ports `simulate_neural_volume.m`.
"""
function simulate_neural_volume(vol_params::VolumeParams,
                                neur_params::NeuronParams,
                                vasc_params::VasculatureParams,
                                dend_params::DendriteParams,
                                axon_params::AxonParams,
                                bg_params::BackgroundParams;
                                rng::AbstractRNG=Random.default_rng())
    finalize!(vol_params)
    vres = vol_params.vres
    vol_sz = vol_params.vol_sz
    H = Int(round(vol_sz[1] * vres))
    W = Int(round(vol_sz[2] * vres))
    D = Int(round(vol_sz[3] * vres))
    vol_depth_vox = Int(round(vol_params.vol_depth * vres))

    # 1. Vasculature.
    if vasc_params.flag != 0
        neur_ves_full, _ = simulate_blood_vessels(vol_params, vasc_params; rng=rng)
    else
        neur_ves_full = falses(H, W, D + vol_depth_vox)
    end

    # 2. Soma positions and meshes.
    locs, Vcells, Vnucs, _ = sample_dense_neurons(neur_params, vol_params,
                                                  neur_ves_full; rng=rng)
    vol_params.N_neur = length(Vcells)

    # 3. Rasterize somata + nuclei.
    neur_soma, neur_vol0, gp_nuc, gp_soma = generate_neural_volume(neur_params,
                                                                    vol_params,
                                                                    locs, Vcells,
                                                                    Vnucs,
                                                                    neur_ves_full)

    # 4. Basal + local-apical dendrites.
    neur_num, dendnum_AD, _ = grow_neuron_dendrites!(vol_params, dend_params,
                                                     neur_soma, neur_ves_full,
                                                     locs, gp_nuc, gp_soma;
                                                     rng=rng)

    # 5. Through-volume apical dendrites.
    neur_num, neur_num_AD = grow_apical_dendrites!(vol_params, dend_params,
                                                   neur_num, dendnum_AD,
                                                   gp_nuc, gp_soma; rng=rng)

    # 6. Per-cell fluorescence.
    gp_vals, neur_vol = set_cell_fluorescence(vol_params, neur_params,
                                              dend_params, neur_num, neur_soma,
                                              neur_num_AD, locs, neur_vol0;
                                              rng=rng)

    # 7. Background neuropil.
    if bg_params.flag != 0
        neur_num, neur_vol, gp_vals, locs =
            generate_bg_dendrites(vol_params, bg_params, dend_params,
                                   neur_vol, neur_num, gp_vals, gp_nuc,
                                   locs; rng=rng)
    end

    # 8. Axons.
    gp_bgvals = Vector{NamedTuple{(:loc, :val),
                                  Tuple{Vector{Int32}, Vector{Float32}}}}()
    bg_proc = Vector{NamedTuple{(:loc, :val),
                                Tuple{Vector{Int32}, Vector{Float32}}}}()
    if axon_params.flag != 0
        neur_vol, gp_bgvals, axon_params = generate_axons(vol_params, axon_params,
                                                          neur_vol, neur_num,
                                                          gp_vals, gp_nuc;
                                                          rng=rng)
        bg_proc = sort_axons(vol_params, axon_params, gp_bgvals,
                             locs .* vol_params.vres; rng=rng)
    end

    # Brain-slab vessel mask (without the above-volume slab).
    neur_ves_brain = falses(H, W, D)
    @inbounds for k in 1:D, j in 1:W, i in 1:H
        neur_ves_brain[i, j, k] = neur_ves_full[i, j, vol_depth_vox + k]
    end

    return NeuralVolume(neur_vol, neur_num, neur_soma, neur_num_AD,
                        gp_nuc, gp_soma, gp_vals, gp_bgvals, bg_proc,
                        neur_ves_brain, neur_ves_full,
                        Matrix{Float64}(locs),
                        collect(Vcells), collect(Vnucs))
end
