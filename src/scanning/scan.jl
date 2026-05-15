# Full-movie scanning with motion + noise.
#
# Ported from upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - scan_volume.m       → scan_volume
#   - imgSubRowShift.m    → img_sub_row_shift
#
# Note: `blurredBackComp2.m` is part of the temporal-focusing scattering
# background and depends on the deferred Chunk-7 cortical-light-path
# orchestrator. Not ported here.

export scan_volume, img_sub_row_shift

"""
    img_sub_row_shift(img_in::AbstractMatrix{<:Real},
                      buf_sz::Integer, x_off, y_off) -> Matrix{Float32}

Sub-select a `(rows - 2·buf_sz, cols - 2·buf_sz)` region of `img_in`
with per-row fractional x and y offsets, modelling fast-axis tissue
motion. `x_off` is either a scalar or a per-row vector; `y_off` is a
per-row vector. Linear interpolation handles fractional shifts; out-of-
bounds rows/columns yield zeros. Ports `imgSubRowShift.m`.
"""
function img_sub_row_shift(img_in::AbstractMatrix{<:Real},
                           buf_sz::Integer, x_off, y_off::AbstractVector)
    H, W = size(img_in)
    nrows = length(y_off)
    x_off_vec = x_off isa AbstractVector ? Float64.(x_off) : fill(Float64(x_off), nrows)
    @assert length(x_off_vec) == nrows
    img_tmp = zeros(Float32, H, W)

    x_off_v = x_off_vec .- buf_sz
    y_off_v = Float64.(y_off) .- buf_sz
    x_off2 = x_off_v .+ (0:(nrows - 1))

    @inbounds for k in 1:nrows
        x2 = x_off2[k]
        if mod(x2, 1) > 0
            lo, hi = floor(Int, x2), ceil(Int, x2)
            t1 = (1 <= lo <= H) ? @view(img_in[lo, :]) : fill(NaN32, W)
            t2 = (1 <= hi <= H) ? @view(img_in[hi, :]) : fill(NaN32, W)
            frac = Float32(mod(x2, 1))
            img_tmp[k, :] .= Float32.(t1) .* (1 - frac) .+ Float32.(t2) .* frac
        else
            xi = Int(x2)
            if 1 <= xi <= H
                img_tmp[k, :] .= Float32.(@view img_in[xi, :])
            else
                img_tmp[k, :] .= NaN32
            end
        end
    end

    offset = max(1, ceil(Int, maximum(abs.(y_off_v))))
    padded = fill(NaN32, H, W + 2 * offset)
    padded[:, (offset + 1):(offset + W)] .= img_tmp

    img_out = zeros(Float32, H, W)
    @inbounds for k in 1:H
        ycol = y_off_v[min(k, end)]
        lo_y, hi_y = floor(Int, ycol), ceil(Int, ycol)
        frac = Float32(mod(ycol, 1))
        for j in 1:W
            t1 = padded[k, lo_y + offset + j]
            t2 = padded[k, hi_y + offset + j]
            img_out[k, j] = isnan(t1) ? Float32(0) :
                            isnan(t2) ? t1 :
                            t1 * (1 - frac) + t2 * frac
        end
    end
    # Trim the buffer rows/cols.
    return img_out[(buf_sz + 1):(H - buf_sz), (buf_sz + 1):(W - buf_sz)]
end

"""
    scan_volume(neur_vol::NeuralVolume, psf::AbstractArray{<:Real,3},
                neur_act, scan_params::ScanParams;
                noise_params=nothing, tpm_params=nothing,
                spike_opts=nothing, rng=Random.default_rng(),
                return_clean::Bool=false, return_motion::Bool=false)
        -> mov  (or (mov, mov_clean) or (mov, mov_clean, mot_hist))

Scan a volume across `nt` time steps to produce a noisy movie. `neur_act`
is a `NamedTuple` of `K × nt` matrices for `soma`/`dend`/`bg`. Returns
the noisy movie; optionally the clean (pre-noise) movie and the per-
frame motion history (`3 × nt`).

Ports `scan_volume.m` minus the temporal-focusing
(`psfT`/`psfB`/`colmask`) and TIFF-streaming output paths.
"""
function scan_volume(neur_vol::NeuralVolume,
                     psf::AbstractArray{<:Real,3},
                     neur_act,
                     scan_params::ScanParams;
                     noise_params=nothing,
                     tpm_params=nothing,
                     spike_opts=nothing,
                     rng::AbstractRNG=Random.default_rng(),
                     return_clean::Bool=false,
                     return_motion::Bool=false)

    H, W, D = size(neur_vol.neur_vol)
    Np1, Np2, Np3 = size(psf)
    sfrac = scan_params.sfrac
    scan_buff = scan_params.scan_buff
    mot_opt = scan_params.motion

    # Activity matrices (K × nt).
    soma_act = Float32.(neur_act.soma)
    dend_act = hasproperty(neur_act, :dend) ? Float32.(neur_act.dend) : copy(soma_act)
    bg_act   = hasproperty(neur_act, :bg)   ? Float32.(neur_act.bg)   :
               zeros(Float32, 1, size(soma_act, 2))
    K_soma = size(soma_act, 1)
    K_axon = size(bg_act, 1)
    nt = size(soma_act, 2)

    sv = setup_scan_volume_frame(neur_vol, psf, scan_params)
    # Use scan_avg + the ScanVolume's cached freq_psf via single_scan path.

    # Sigscale: TPM signal scaling × dt × sfrac² / volume xy area.
    sigscale = if !isnothing(tpm_params) && !isnothing(spike_opts)
        Float32(tpm_signal_scale(tpm_params) * spike_opts.dt * sfrac^2 /
                (250 * 250))
    else
        1.0f0
    end

    # Output movie size after binning + buffer crop.
    H_out = (H - 2 * scan_buff) ÷ sfrac
    W_out = (W - 2 * scan_buff) ÷ sfrac
    mov      = zeros(Float32, H_out, W_out, nt)
    mov_clean = return_clean ? zeros(Float32, H_out, W_out, nt) : Array{Float32}(undef, 0, 0, 0)
    mot_hist = return_motion ? zeros(Float32, 3, nt) : Array{Float32}(undef, 0, 0)

    # Motion parameters (matching upstream).
    z_base = Int(floor(0.5 * (D - Np3)))
    x_loc = scan_buff + 1
    y_loc = scan_buff + 1
    z_loc = z_base

    if mot_opt
        d_stps   = vcat([-1, 1], zeros(Int, 5))
        d_stpsZ  = vcat([-1, 1], zeros(Int, 100))
        d_stps2  = collect(-3:3)
        p_jump   = 0.05
        maxshear = 1.0 / 200
        zmaxdiff = 2
    else
        d_stps   = [0]
        d_stpsZ  = [0]
        d_stps2  = [0]
        p_jump   = 0.0
        maxshear = 0.0
        zmaxdiff = 0
    end

    # Per-cell soma/dend/axon voxel/value vectors (already split in sv).
    cutoff = 1e-2

    for kk in 1:nt
        if rand(rng) > p_jump
            x_loc = clamp(x_loc + rand(rng, d_stps2), 1, 2 * scan_buff + 1)
            y_loc = clamp(y_loc + rand(rng, d_stps2), 1, 2 * scan_buff + 1)
        end
        x_pos = clamp(x_loc + rand(rng, d_stps), 1, 2 * scan_buff + 1)
        y_pos = clamp(y_loc + rand(rng, d_stps), 1, 2 * scan_buff + 1)
        z_loc = clamp(z_loc + rand(rng, d_stpsZ), z_base - zmaxdiff,
                      z_base + zmaxdiff)
        z_loc = clamp(z_loc + rand(rng, d_stps), 1, D - Np3 + 1)
        return_motion && (mot_hist[:, kk] .= Float32[x_pos, y_pos, z_loc])

        # Per-row y shear vector.
        y_shr = zeros(Float32, H)
        if mot_opt
            head_len = rand(rng, 1:max(1, fld(2 * H, 5)))
            mid_len  = max(1, Int(round(rand(rng) * 3 * H / 5)))
            slope = (2 * (rand(rng) - 0.5)) * maxshear * H
            seg = collect(range(0, 1; length=mid_len)) .* slope
            for i in 1:H
                if i <= head_len
                    y_shr[i] = 0
                elseif i - head_len <= mid_len
                    y_shr[i] = Float32(seg[i - head_len])
                else
                    y_shr[i] = Float32(seg[end])
                end
            end
        end
        y_off = clamp.(y_pos .+ y_shr .+ rand(rng, d_stps, H), 1, 2 * scan_buff + 1)

        # Build TMPvol with this frame's activity.
        TMPvol = zeros(Float32, H, W, D)
        @inbounds for ll in 1:K_soma
            a = soma_act[ll, kk]
            a > cutoff || continue
            isempty(sv.soma_loc[ll]) && continue
            for (i, li) in enumerate(sv.soma_loc[ll])
                TMPvol[Int(li)] = sv.soma_val[ll][i] * a
            end
        end
        @inbounds for ll in 1:K_soma
            a = dend_act[ll, kk]
            a > cutoff || continue
            isempty(sv.dend_loc[ll]) && continue
            for (i, li) in enumerate(sv.dend_loc[ll])
                TMPvol[Int(li)] = sv.dend_val[ll][i] * a
            end
        end
        @inbounds for ll in 1:K_axon
            a = bg_act[ll, kk]
            a > cutoff || continue
            ll <= length(sv.axon_loc) || continue
            isempty(sv.axon_loc[ll]) && continue
            for (i, li) in enumerate(sv.axon_loc[ll])
                TMPvol[Int(li)] += sv.axon_val[ll][i] * a
            end
        end

        # Scan the active slab.
        z_hi = min(z_loc + Np3 - 1, D)
        sub_vol = @view TMPvol[:, :, z_loc:z_hi]
        clean_img = (sigscale / (2 * sfrac^2)) .*
                    single_scan(sub_vol, (Np1, Np2, Np3), sv.freq_psf;
                                z_sub=scan_params.scan_avg, freq_opt=true)

        # Per-row shift to model motion + buffer crop.
        if mot_opt && size(clean_img, 1) > 2 * scan_buff && size(clean_img, 2) > 2 * scan_buff
            clean_img = img_sub_row_shift(clean_img, scan_buff, x_pos, round.(Int, y_off))
        else
            clean_img = clean_img[(scan_buff + 1):(size(clean_img, 1) - scan_buff),
                                   (scan_buff + 1):(size(clean_img, 2) - scan_buff)]
        end

        if sfrac > 1 && floor(sfrac) == sfrac
            sf = Int(sfrac)
            H2 = size(clean_img, 1) ÷ sf
            W2 = size(clean_img, 2) ÷ sf
            binned = zeros(Float32, H2, W2)
            for j in 1:W2, i in 1:H2
                s = 0f0
                @inbounds for jj in 1:sf, ii in 1:sf
                    s += clean_img[(i - 1) * sf + ii, (j - 1) * sf + jj]
                end
                binned[i, j] = s
            end
            clean_img = binned
        end

        if return_clean
            # Trim to expected output size.
            h_trim = min(size(clean_img, 1), H_out)
            w_trim = min(size(clean_img, 2), W_out)
            mov_clean[1:h_trim, 1:w_trim, kk] .= @view clean_img[1:h_trim, 1:w_trim]
        end

        # Noise model.
        samp_img = if !isnothing(noise_params)
            np = noise_params.sigscale != noise_params.sigscale ? noise_params :
                noise_params  # no-op; sigscale already encoded above
            noisy = poisson_gauss_noise(clean_img, noise_params; rng=rng)
            pixel_bleed(noisy, noise_params.bleedp, noise_params.bleedw; rng=rng)
        else
            clean_img
        end

        h_trim = min(size(samp_img, 1), H_out)
        w_trim = min(size(samp_img, 2), W_out)
        mov[1:h_trim, 1:w_trim, kk] .= @view samp_img[1:h_trim, 1:w_trim]
    end

    if return_clean && return_motion
        return mov, mov_clean, mot_hist
    elseif return_clean
        return mov, mov_clean
    elseif return_motion
        return mov, mot_hist
    else
        return mov
    end
end
