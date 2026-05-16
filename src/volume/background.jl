# Background-dendrite (neuropil) generation.
#
# Ported from the upstream NAOMi-Sim file (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - generate_bgdendrites.m → generate_bg_dendrites

export generate_bg_dendrites

"""
    generate_bg_dendrites(vol_params, bg_params, dend_params, neur_vol,
                          neur_num, gp_vals, gp_nuc, neur_locs;
                          neur_vol_flag=true, rng=Random.default_rng())
        -> (neur_num_out, neur_vol_out, gp_vals_out, neur_locs_out)

Generate background-dendrite processes rooted *outside* the neural
volume and walking *into* it via [`dendrite_random_walk`](@ref). The
resulting voxels are dilated via [`dilate_dendrite_paths_all`](@ref)
and merged into `neur_num` with labels starting at
`vol_params.N_neur + vol_params.N_den + 1`. Per-voxel fluorescence
follows the dendrite radial-decay rule
`(wtSc[2]·exp(−1) + (1 − wtSc[2]))·(1 − wtSc[3]·rand)`.

Returns the updated `neur_num`, `neur_vol`, `gp_vals` (with new
entries appended for each new background process), and `neur_locs`
(with the root of each new process appended in micron units).

Ports `generate_bgdendrites.m`.
"""
function generate_bg_dendrites(vol_params::VolumeParams,
                               bg_params::BackgroundParams,
                               dend_params::DendriteParams,
                               neur_vol::AbstractArray{<:Real,3},
                               neur_num::AbstractArray{<:Integer,3},
                               gp_vals::AbstractVector,
                               gp_nuc::AbstractVector,
                               neur_locs::AbstractMatrix{<:Real};
                               neur_vol_flag::Bool=true,
                               rng::AbstractRNG=Random.default_rng())
    vres = vol_params.vres
    vol_sz = vol_params.vol_sz
    H = Int(round(vol_sz[1] * vres))
    W = Int(round(vol_sz[2] * vres))
    D = Int(round(vol_sz[3] * vres))
    volsize = (H, W, D)

    bg_pix = falses(volsize)
    for li in eachindex(neur_num)
        bg_pix[li] = neur_num[li] == 0
    end
    for kk in 1:vol_params.N_neur
        for x in gp_nuc[kk][1]
            bg_pix[Int(x)] = false
        end
    end

    dtParams = collect(Float64, dend_params.dtParams)
    dtParams[2:3] .*= vres
    thicknessScale = dend_params.thicknessScale * vres * vres

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

    M = rand(rng, Float32, volsize)
    for li in eachindex(M)
        if !bg_pix[li]
            M[li] = floatmax(Float32)
        end
    end
    M[1, :, :]   .= floatmax(Float32)
    M[end, :, :] .= floatmax(Float32)
    M[:, 1, :]   .= floatmax(Float32)
    M[:, end, :] .= floatmax(Float32)
    M[:, :, 1]   .= floatmax(Float32)
    M[:, :, end] .= floatmax(Float32)

    dendVar = 0.25
    dtSize = (dtParams[2], dtParams[2], dtParams[3])
    idxvol = zeros(UInt16, volsize)
    numvol = zeros(Float32, volsize)

    maxlength  = bg_params.maxlength
    distsc     = bg_params.distsc
    fillweight = bg_params.fillweight
    maxel      = bg_params.maxel
    minlength  = bg_params.minlength
    shiftdist  = 3

    # Number of processes to generate (matches upstream's product-ratio rule).
    Nproc_target = Int(round(((prod(volsize .+ 2 .* dtSize) / prod(volsize)) - 1) *
                              vol_params.N_neur))
    idx = 0
    neur_locs_out = collect(Float64, neur_locs)
    new_locs = Float64[]

    for j in 1:Nproc_target
        # Root must be outside the volume in at least one axis.
        root = (Int(floor(rand(rng) * (volsize[1] + 2 * dtSize[1]) - dtSize[1])),
                Int(floor(rand(rng) * (volsize[2] + 2 * dtSize[2]) - dtSize[2])),
                Int(floor(rand(rng) * (volsize[3] + 2 * dtSize[3]) - dtSize[3])))
        tries = 0
        while (1 <= root[1] <= volsize[1] && 1 <= root[2] <= volsize[2] &&
               1 <= root[3] <= volsize[3]) && tries < 100
            tries += 1
            root = (Int(floor(rand(rng) * (volsize[1] + 2 * dtSize[1]) - dtSize[1])),
                    Int(floor(rand(rng) * (volsize[2] + 2 * dtSize[2]) - dtSize[2])),
                    Int(floor(rand(rng) * (volsize[3] + 2 * dtSize[3]) - dtSize[3])))
        end
        # Record root micron position.
        push!(new_locs, root[1] / vres)
        push!(new_locs, root[2] / vres)
        push!(new_locs, root[3] / vres)

        dendpts_lin = Int32[]
        for _ in 1:Int(round(dtParams[1]))
            θ = rand(rng) * 2π
            r = sqrt(rand(rng)) * dtParams[2]
            dends = (Int(floor(r * cos(θ) + root[1])),
                     Int(floor(r * sin(θ) + root[2])),
                     Int(floor(2 * dtParams[3] * (rand(rng) - 0.5) + root[3])))
            (1 <= dends[1] <= volsize[1] && 1 <= dends[2] <= volsize[2] &&
             1 <= dends[3] <= volsize[3]) || continue
            # Find an entry point on the volume boundary along (root → dends).
            v_root = collect(Float64, root)
            v_dends = collect(Float64, dends)
            diff = v_dends .- v_root
            best_alpha = 0.0
            shiftLoc = 1
            for ax in 1:3
                if v_root[ax] < 1
                    α = (1 - v_root[ax]) / diff[ax]
                    if α > best_alpha
                        best_alpha = α
                        shiftLoc = ax
                    end
                end
                if v_root[ax] > volsize[ax]
                    α = (volsize[ax] - v_root[ax]) / diff[ax]
                    if α > best_alpha
                        best_alpha = α
                        shiftLoc = ax + 3
                    end
                end
            end
            root2 = round.(Int, best_alpha .* diff .+ v_root)
            # Jitter perpendicular to the entry axis.
            j1 = max(1, ceil(Int, shiftdist * rand(rng)))
            j2 = max(1, ceil(Int, shiftdist * rand(rng)))
            if shiftLoc == 1 || shiftLoc == 4
                root2[2] += j1; root2[3] += j2
            elseif shiftLoc == 2 || shiftLoc == 5
                root2[1] += j1; root2[3] += j2
            else
                root2[1] += j1; root2[2] += j2
            end
            for ax in 1:3
                root2[ax] = clamp(root2[ax], 1, volsize[ax])
            end
            bgpts = NTuple{3,Int}[]
            numit = 0
            while isempty(bgpts) && numit < 30
                numit += 1
                bgpts = dendrite_random_walk(M, (root2[1], root2[2], root2[3]), dends;
                                             distsc=distsc,
                                             maxlength=Int(maxlength),
                                             fillweight=fillweight,
                                             maxel=maxel,
                                             minlength=Int(minlength))
            end
            if !isempty(bgpts)
                pushfirst!(bgpts, (root2[1], root2[2], root2[3]))
                dendSz = max(0.0, 1 + dendVar * randn(rng))^2
                n = length(bgpts)
                pw = ones(Float32, n)
                if n > 2
                    for k in 2:(n - 1)
                        p, c, q = bgpts[k - 1], bgpts[k], bgpts[k + 1]
                        d2 = abs(2c[1] - p[1] - q[1]) +
                             abs(2c[2] - p[2] - q[2]) +
                             abs(2c[3] - p[3] - q[3])
                        pw[k] = Float32(dendSz * (1 - (1 - 1 / sqrt(2)) * d2 / 2))
                    end
                    pw[1] = pw[2]; pw[end] = pw[end - 1]
                    pw .= max.(pw, 0f0)
                else
                    pw .= Float32(dendSz)
                end
                lin = LinearIndices(volsize)
                for (k, p) in enumerate(bgpts)
                    li = Int32(lin[p...])
                    push!(dendpts_lin, li)
                    numvol[Int(li)] = pw[k]
                end
            end
        end

        if !isempty(dendpts_lin)
            idx += 1
            for li in dendpts_lin
                idxvol[Int(li)] = UInt16(idx)
                numvol[Int(li)] *= Float32(thicknessScale * dtParams[4])
            end
        end
    end

    _, pathnum = dilate_dendrite_paths_all(numvol, idxvol, .!bg_pix; rng=rng)
    Ncomps = vol_params.N_neur + Int(round(vol_params.N_den))
    pathnum_offset = copy(pathnum)
    for li in eachindex(pathnum_offset)
        if pathnum_offset[li] > 0
            pathnum_offset[li] += Ncomps
        end
    end
    neur_num_out = UInt16.(neur_num) .+ UInt16.(pathnum_offset)

    wtSc = collect(Float64, dend_params.weightScale)
    gp_vals_out = copy(gp_vals)
    for i in (Ncomps + 1):(Ncomps + idx)
        cidx = Int32[Int32(li) for li in eachindex(neur_num_out)
                     if neur_num_out[li] == UInt16(i)]
        vals = Float32[(wtSc[2] * exp(-(dtParams[2] / vres) / wtSc[1]) +
                        (1 - wtSc[2])) * (1 - wtSc[3] * rand(rng))
                       for _ in 1:length(cidx)]
        push!(gp_vals_out, (loc=cidx, val=vals, is_soma=BitVector(falses(length(cidx)))))
        if neur_vol_flag
            for (k, li) in enumerate(cidx)
                neur_vol_out[Int(li)] = vals[k]
            end
        end
    end

    # Append generated background-process roots to neur_locs (micron units).
    if !isempty(new_locs)
        n_new = length(new_locs) ÷ 3
        appended = reshape(new_locs, 3, n_new)'
        neur_locs_out = vcat(neur_locs_out, appended[1:min(idx, n_new), :])
    end

    return neur_num_out, neur_vol_out, gp_vals_out, neur_locs_out
end
