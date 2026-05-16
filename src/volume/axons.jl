# Axon / neuropil-background generation.
#
# Ported from the upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - generate_axons.m → generate_axons
#   - sort_axons.m     → sort_axons

export generate_axons, sort_axons

"""
    generate_axons(vol_params, axon_params, neur_vol, neur_num, gp_vals, gp_nuc;
                   neur_vol_flag=true, rng=Random.default_rng())
        -> (neur_vol_out, gp_bgvals, axon_params_out)

Fill background voxels with many short directed-random-walk processes
("axons"). Each process has a small random number of branches; per-voxel
values are `1/axon_params.maxel · (1 + axon_params.varfill·N(0,1))`.

Returns:
- `neur_vol_out::Array{Float32,3}` — fluorescence map with axons added
  (unchanged if `neur_vol_flag = false`).
- `gp_bgvals::Vector{NamedTuple{(:loc, :val)}}` — per-process voxel
  index lists and values.
- `axon_params_out` — possibly updated `axon_params` (preserves a
  modified `N_proc` for downstream `sort_axons`).

Ports `generate_axons.m`.
"""
function generate_axons(vol_params::VolumeParams,
                        axon_params::AxonParams,
                        neur_vol::AbstractArray{<:Real,3},
                        neur_num::AbstractArray{<:Integer,3},
                        gp_vals::AbstractVector,
                        gp_nuc::AbstractVector;
                        neur_vol_flag::Bool=true,
                        rng::AbstractRNG=Random.default_rng())
    vol_sz = vol_params.vol_sz
    vres = vol_params.vres
    H = Int(round(vol_sz[1] * vres))
    W = Int(round(vol_sz[2] * vres))
    D = Int(round(vol_sz[3] * vres))
    volsize = (H, W, D)

    bg_pix = falses(volsize)
    for li in eachindex(neur_num)
        bg_pix[li] = neur_num[li] == 0
    end
    for kk in 1:length(gp_nuc)
        for x in gp_nuc[kk][1]
            bg_pix[Int(x)] = false
        end
    end

    fillnum = round(Int, axon_params.maxfill * axon_params.maxvoxel * sum(bg_pix))
    N_bg = vol_params.N_bg

    neur_vol_out = if neur_vol_flag
        v = zeros(Float32, volsize)
        for kk in 1:length(gp_vals)
            for (i, li) in enumerate(gp_vals[kk].loc)
                v[Int(li)] = gp_vals[kk].val[i]
            end
            if kk <= length(gp_nuc)
                idxs, val = gp_nuc[kk]
                for x in idxs
                    v[Int(x)] = Float32(val)
                end
            end
        end
        v
    else
        Float32.(neur_vol)
    end

    padsize = axon_params.padsize
    volpad = (H + 2 * padsize, W + 2 * padsize, D + 2 * padsize)
    M = rand(rng, Float32, volpad)
    # Mark non-background voxels as blocked (max float).
    for k in 1:D, j in 1:W, i in 1:H
        if !bg_pix[i, j, k]
            M[i + padsize, j + padsize, k + padsize] = floatmax(Float32)
        end
    end

    gp_bgvals = NamedTuple{(:loc, :val), Tuple{Vector{Int32}, Vector{Float32}}}[]
    j = 0
    numit2 = 0
    nummax = 10_000
    maxfill_val = axon_params.fillweight * axon_params.maxvoxel
    threshold_max = floatmax(Float32) / 2

    while fillnum > 0 && j < N_bg && numit2 < nummax
        bgpts = NTuple{3,Int}[]
        numit2 = 0
        while length(bgpts) < axon_params.minlength && numit2 < nummax
            numit2 += 1
            root = (max(1, ceil(Int, (volpad[1] - 2) * rand(rng) + 1)),
                    max(1, ceil(Int, (volpad[2] - 2) * rand(rng) + 1)),
                    max(1, ceil(Int, (volpad[3] - 2) * rand(rng) + 1)))
            tries = 0
            while M[root...] > maxfill_val && tries < 100
                tries += 1
                root = (max(1, ceil(Int, (volpad[1] - 2) * rand(rng) + 1)),
                        max(1, ceil(Int, (volpad[2] - 2) * rand(rng) + 1)),
                        max(1, ceil(Int, (volpad[3] - 2) * rand(rng) + 1)))
            end
            tries == 100 && break
            ends = (clamp(Int(ceil(root[1] + 2 * axon_params.maxdist * vres *
                                       (rand(rng) - 0.5))), 1, volpad[1]),
                    clamp(Int(ceil(root[2] + 2 * axon_params.maxdist * vres *
                                       (rand(rng) - 0.5))), 1, volpad[2]),
                    clamp(Int(ceil(root[3] + 2 * axon_params.maxdist * vres *
                                       (rand(rng) - 0.5))), 1, volpad[3]))
            bgpts = dendrite_random_walk(M, root, ends;
                                         distsc=axon_params.distsc,
                                         maxlength=Int(axon_params.maxlength),
                                         fillweight=axon_params.fillweight,
                                         maxel=axon_params.maxvoxel,
                                         minlength=Int(axon_params.minlength))
        end
        numit2 >= nummax && break

        # Add branches off the main trunk.
        nbranches = max(0, Int(round(axon_params.numbranches +
                                     axon_params.varbranches * randn(rng))))
        for _ in 1:nbranches
            bgpts2 = NTuple{3,Int}[]
            numit = 0
            while length(bgpts2) < axon_params.minlength && numit < 100
                numit += 1
                root = bgpts[max(1, ceil(Int, length(bgpts) * rand(rng)))]
                # avoid roots on the volume boundary
                tries = 0
                while (root[1] == 1 || root[1] == volpad[1] ||
                       root[2] == 1 || root[2] == volpad[2] ||
                       root[3] == 1 || root[3] == volpad[3]) && tries < 100
                    tries += 1
                    root = bgpts[max(1, ceil(Int, length(bgpts) * rand(rng)))]
                end
                ends = (clamp(Int(ceil(root[1] + 2 * axon_params.maxdist * vres *
                                       (rand(rng) - 0.5))), 1, volpad[1]),
                        clamp(Int(ceil(root[2] + 2 * axon_params.maxdist * vres *
                                       (rand(rng) - 0.5))), 1, volpad[2]),
                        clamp(Int(ceil(root[3] + 2 * axon_params.maxdist * vres *
                                       (rand(rng) - 0.5))), 1, volpad[3]))
                bgpts2 = dendrite_random_walk(M, root, ends;
                                              distsc=axon_params.distsc,
                                              maxlength=Int(axon_params.maxlength),
                                              fillweight=axon_params.fillweight,
                                              maxel=axon_params.maxvoxel,
                                              minlength=Int(axon_params.minlength))
            end
            append!(bgpts, bgpts2)
        end

        # Translate back from padded coordinates.
        translated = NTuple{3,Int}[]
        for p in bgpts
            q = (p[1] - padsize, p[2] - padsize, p[3] - padsize)
            (1 <= q[1] <= H && 1 <= q[2] <= W && 1 <= q[3] <= D) || continue
            push!(translated, q)
        end

        if !isempty(translated)
            lin = LinearIndices(volsize)
            loc = Int32[Int32(lin[p...]) for p in translated]
            vals = fill(Float32((1 / axon_params.maxel) *
                                max(0, 1 + axon_params.varfill * randn(rng))),
                        length(translated))
            push!(gp_bgvals, (loc=loc, val=vals))
            fillnum -= length(translated)
            j += 1
            if neur_vol_flag
                for (i, li) in enumerate(loc)
                    neur_vol_out[Int(li)] += vals[i]
                end
            end
        end
    end

    new_params = AxonParams(flag=axon_params.flag,
                            distsc=axon_params.distsc,
                            fillweight=axon_params.fillweight,
                            maxlength=axon_params.maxlength,
                            minlength=axon_params.minlength,
                            maxdist=axon_params.maxdist,
                            maxel=axon_params.maxel,
                            varfill=axon_params.varfill,
                            maxvoxel=axon_params.maxvoxel,
                            padsize=axon_params.padsize,
                            numbranches=axon_params.numbranches,
                            varbranches=axon_params.varbranches,
                            maxfill=axon_params.maxfill,
                            N_proc=j,
                            l=axon_params.l,
                            rho=axon_params.rho)
    return neur_vol_out, gp_bgvals, new_params
end

"""
    sort_axons(vol_params, axon_params, gp_bgvals, cell_pos;
               rng=Random.default_rng())
        -> Vector{NamedTuple{(:loc, :val)}}

Bin the axon processes returned by [`generate_axons`](@ref) into
`axon_params.N_proc` groups. If `N_proc > N_neur + N_den`, each cell
gets paired with its spatially-nearest axon process (greedy
assignment); the remainder are dropped into random extra bins.
Otherwise axons are binned uniformly at random.

Returns one entry per group, each with concatenated `(loc, val)`
vectors. Ports `sort_axons.m`.
"""
function sort_axons(vol_params::VolumeParams,
                    axon_params::AxonParams,
                    gp_bgvals::AbstractVector,
                    cell_pos::AbstractMatrix{<:Real};
                    rng::AbstractRNG=Random.default_rng())
    N_proc = axon_params.N_proc
    vol_sz = vol_params.vol_sz
    vres = vol_params.vres
    H = Int(round(vol_sz[1] * vres))
    W = Int(round(vol_sz[2] * vres))
    D = Int(round(vol_sz[3] * vres))
    cart = CartesianIndices((H, W, D))

    bg_proc = [(loc=Int32[], val=Float32[]) for _ in 1:N_proc]
    # accumulate via mutable wrappers
    accum = [(loc=Int32[], val=Float32[]) for _ in 1:N_proc]

    if N_proc > vol_params.N_neur + Int(round(vol_params.N_den))
        N_comps = vol_params.N_neur + Int(round(vol_params.N_den))
        nbg = length(gp_bgvals)
        gp_bgpos = zeros(Float64, nbg, 3)
        for kk in 1:nbg
            e = gp_bgvals[kk]
            isempty(e.loc) && continue
            xs = 0.0; ys = 0.0; zs = 0.0
            for li in e.loc
                c = cart[Int(li)]
                xs += c[1]; ys += c[2]; zs += c[3]
            end
            n = length(e.loc)
            gp_bgpos[kk, :] = [xs / n, ys / n, zs / n]
        end
        # Distance matrix between first N_comps cells and bg processes.
        cp = view(cell_pos, 1:min(N_comps, size(cell_pos, 1)), :)
        dist_mat = fill(Inf, size(cp, 1), nbg)
        for i in 1:size(cp, 1), j in 1:nbg
            dist_mat[i, j] = sqrt((cp[i, 1] - gp_bgpos[j, 1])^2 +
                                  (cp[i, 2] - gp_bgpos[j, 2])^2 +
                                  (cp[i, 3] - gp_bgpos[j, 3])^2)
        end
        idxlist = falses(nbg)
        for ii in 1:size(cp, 1)
            (_, idx) = findmin(@view dist_mat[ii, :])
            dist_mat[:, idx] .= Inf
            accum[ii] = gp_bgvals[idx]
            idxlist[idx] = true
        end
        for kk in 1:nbg
            idxlist[kk] && continue
            extra = N_comps + max(1, ceil(Int, (N_proc - N_comps) * rand(rng)))
            cur = accum[extra]
            accum[extra] = (loc=vcat(cur.loc, gp_bgvals[kk].loc),
                            val=vcat(cur.val, gp_bgvals[kk].val))
        end
    else
        for kk in 1:length(gp_bgvals)
            idx = max(1, ceil(Int, N_proc * rand(rng)))
            cur = accum[idx]
            accum[idx] = (loc=vcat(cur.loc, gp_bgvals[kk].loc),
                          val=vcat(cur.val, gp_bgvals[kk].val))
        end
    end
    return accum
end
