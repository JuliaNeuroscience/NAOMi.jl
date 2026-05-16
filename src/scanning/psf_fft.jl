# PSF FFT + single-frame scan.
#
# Ported from upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - psf_fft.m         → psf_fft
#   - single_scan.m     → single_scan
#   - setup_scan_volume_frame.m → setup_scan_volume_frame
#   - scan_volume_frame.m       → scan_volume_frame

export psf_fft, single_scan, setup_scan_volume_frame, scan_volume_frame,
       ScanVolume

"""
    psf_fft(vol_sz::NTuple{3,<:Integer}, psf::AbstractArray{<:Real,3};
            z_sub::Integer=1) -> Array{ComplexF32,3}

Pre-compute the 2-D FFT (along x and y) of `psf` zero-padded to
`vol_sz + size(psf) - 1` in the x/y plane, optionally pre-summing every
`z_sub` axial slices. The result is fed to [`single_scan`](@ref) as
`freq_psf`. Ports `psf_fft.m`.
"""
function psf_fft(vol_sz::NTuple{3,<:Integer}, psf::AbstractArray{<:Real,3};
                 z_sub::Integer=1)
    if z_sub > 1
        N3 = size(psf, 3)
        N_slce = cld(N3, z_sub)
        psf2 = zeros(Float32, size(psf, 1), size(psf, 2), N_slce)
        for k in 1:N_slce
            base = (k - 1) * z_sub + 1
            for off in 0:(z_sub - 1)
                slc = base + off
                if slc <= N3
                    psf2[:, :, k] .+= Float32.(@view psf[:, :, slc])
                end
            end
        end
        sz = (vol_sz[1] + size(psf2, 1) - 1, vol_sz[2] + size(psf2, 2) - 1)
        padded = zeros(ComplexF32, sz[1], sz[2], size(psf2, 3))
        for k in axes(psf2, 3)
            padded[1:size(psf2, 1), 1:size(psf2, 2), k] .= ComplexF32.(@view psf2[:, :, k])
        end
        return fft(padded, (1, 2))
    else
        sz = (vol_sz[1] + size(psf, 1) - 1, vol_sz[2] + size(psf, 2) - 1)
        padded = zeros(ComplexF32, sz[1], sz[2], size(psf, 3))
        for k in axes(psf, 3)
            padded[1:size(psf, 1), 1:size(psf, 2), k] .= ComplexF32.(@view psf[:, :, k])
        end
        return fft(padded, (1, 2))
    end
end

"""
    single_scan(neur_vol::AbstractArray{<:Real,3}, psf_sz::NTuple{3,<:Integer},
                psf_or_freq::AbstractArray; z_sub::Integer=1,
                freq_opt::Bool=false, fwd_plan=nothing, inv_plan=nothing)
        -> Matrix{Float32}

Convolve a 3-D fluorescence volume with the supplied PSF (either spatial
or pre-FFTed). Returns a 2-D image cropped to the volume's xy extent.
Ports `single_scan.m`. With `freq_opt=true`, `psf_or_freq` must be the
output of [`psf_fft`](@ref).

`fwd_plan` / `inv_plan` are optional pre-built FFT plans (see
`AbstractFFTs.plan_fft`) for the frequency-domain path: `fwd_plan` for a
batched 2-D transform of a `size(psf_or_freq)` array, `inv_plan` for the
inverse 2-D transform of its `(:, :)` slice. Supplying them lets a caller
that scans many frames avoid re-planning on every call.
"""
function single_scan(neur_vol::AbstractArray{<:Real,3},
                     psf_sz::NTuple{3,<:Integer},
                     psf_or_freq::AbstractArray;
                     z_sub::Integer=1, freq_opt::Bool=false,
                     fwd_plan=nothing, inv_plan=nothing)
    N3 = size(neur_vol, 3)
    if z_sub > 1
        N_slce = cld(N3, z_sub)
        nv2 = zeros(Float32, size(neur_vol, 1), size(neur_vol, 2), N_slce)
        for k in 1:N_slce
            base = (k - 1) * z_sub + 1
            for off in 0:(z_sub - 1)
                slc = base + off
                if slc <= N3
                    nv2[:, :, k] .+= Float32.(@view neur_vol[:, :, slc])
                end
            end
        end
        if !freq_opt
            psf2 = zeros(Float32, size(psf_or_freq, 1), size(psf_or_freq, 2), N_slce)
            for k in 1:N_slce
                base = (k - 1) * z_sub + 1
                for off in 0:(z_sub - 1)
                    slc = base + off
                    if slc <= size(psf_or_freq, 3)
                        psf2[:, :, k] .+= Float32.(@view psf_or_freq[:, :, slc])
                    end
                end
            end
        end
        neur_vol_eff = nv2
    else
        neur_vol_eff = neur_vol
        if !freq_opt
            psf2 = psf_or_freq
        end
    end

    if freq_opt
        sz1, sz2, Npz = size(psf_or_freq)
        H1, W1 = size(neur_vol_eff, 1), size(neur_vol_eff, 2)
        nz = min(size(neur_vol_eff, 3), Npz)
        # Pad every axial slice into one (sz1, sz2, Npz) array and run a
        # single batched 2-D FFT over all slices, then convolve with the
        # pre-FFTed PSF and sum along z.
        slabs = zeros(ComplexF32, sz1, sz2, Npz)
        for k in 1:nz
            @view(slabs[1:H1, 1:W1, k]) .= ComplexF32.(@view neur_vol_eff[:, :, k])
        end
        fS = fwd_plan === nothing ? fft(slabs, (1, 2)) : fwd_plan * slabs
        acc = zeros(ComplexF32, sz1, sz2)
        for k in 1:nz
            @views acc .+= fS[:, :, k] .* psf_or_freq[:, :, k]
        end
        img_full = real.(inv_plan === nothing ? ifft(acc, (1, 2)) : inv_plan * acc)
        y_ix = cld(psf_sz[1] - 1, 2) .+ (1, H1)
        y_jx = cld(psf_sz[2] - 1, 2) .+ (1, W1)
        return Float32.(img_full[y_ix[1]:y_ix[2], y_jx[1]:y_jx[2]])
    else
        # Spatial-domain conv2 per slice, "same" cropping.
        img = zeros(Float32, size(neur_vol_eff, 1), size(neur_vol_eff, 2))
        psf3 = freq_opt ? nothing : psf2
        nz = min(size(neur_vol_eff, 3), size(psf3, 3))
        for k in 1:nz
            img .+= _conv2_same(@view(neur_vol_eff[:, :, k]),
                                @view(psf3[:, :, k]))
        end
        return img
    end
end

# 2-D "same"-mode convolution (zero-padded boundaries). Hand-rolled to
# match MATLAB `conv2(A, B, 'same')`.
function _conv2_same(A::AbstractMatrix, B::AbstractMatrix)
    H, W = size(A)
    h, w = size(B)
    pad_h = cld(h - 1, 2)
    pad_w = cld(w - 1, 2)
    out = zeros(Float32, H, W)
    for j in 1:W, i in 1:H
        s = 0.0f0
        for jj in 1:w, ii in 1:h
            ai = i + ii - 1 - pad_h
            aj = j + jj - 1 - pad_w
            (1 <= ai <= H && 1 <= aj <= W) || continue
            s += Float32(A[ai, aj]) * Float32(B[h - ii + 1, w - jj + 1])
        end
        out[i, j] = s
    end
    return out
end

"""
    ScanVolume

Pre-processed per-cell soma/dendrite/axon voxel-and-value arrays used by
[`scan_volume_frame`](@ref) to assemble a fluorescence volume per frame.

Built by [`setup_scan_volume_frame`](@ref).
"""
struct ScanVolume
    soma_loc::Vector{Vector{Int32}}
    soma_val::Vector{Vector{Float32}}
    dend_loc::Vector{Vector{Int32}}
    dend_val::Vector{Vector{Float32}}
    axon_loc::Vector{Vector{Int32}}
    axon_val::Vector{Vector{Float32}}
    nuc_loc::Vector{Vector{Int32}}
    nuc_val::Vector{Vector{Float32}}
    freq_psf::Array{ComplexF32, 3}
    psf_sz::NTuple{3, Int}
    vol_sz::NTuple{3, Int}
    g_blur::Union{Nothing, Matrix{Float32}}
    nuc_label::Bool
end

"""
    setup_scan_volume_frame(neur_vol::NeuralVolume, psf::AbstractArray{<:Real,3},
                            scan_params::ScanParams;
                            nuc_label::Bool=false, g_blur=nothing) -> ScanVolume

Build the pre-processed [`ScanVolume`](@ref) struct needed for per-frame
scanning. Ports `setup_scan_volume_frame.m` *without* temporal-focusing
psfT/psfB / colmask paths (these are deferred to the cortical-light-path
orchestrator port).
"""
function setup_scan_volume_frame(neur_vol::NeuralVolume,
                                 psf::AbstractArray{<:Real,3},
                                 scan_params::ScanParams;
                                 nuc_label::Bool=false,
                                 g_blur=nothing)
    H, W, D = size(neur_vol.neur_vol)
    Np1, Np2, Np3 = size(psf)
    (H >= Np1 && W >= Np2) || throw(ArgumentError("PSF lateral extent exceeds volume"))
    D >= Np3 || throw(ArgumentError("PSF depth exceeds volume depth"))

    # Split gp_vals into per-cell soma and dendrite vectors.
    K = length(neur_vol.gp_vals)
    soma_loc = [Int32[] for _ in 1:K]
    soma_val = [Float32[] for _ in 1:K]
    dend_loc = [Int32[] for _ in 1:K]
    dend_val = [Float32[] for _ in 1:K]
    for kk in 1:K
        e = neur_vol.gp_vals[kk]
        for (i, li) in enumerate(e.loc)
            if e.is_soma[i]
                push!(soma_loc[kk], li)
                push!(soma_val[kk], e.val[i])
            else
                push!(dend_loc[kk], li)
                push!(dend_val[kk], e.val[i])
            end
        end
    end

    # Axon (background) processes.
    axon_loc = Vector{Int32}[]
    axon_val = Vector{Float32}[]
    for e in neur_vol.bg_proc
        push!(axon_loc, Int32.(e.loc))
        push!(axon_val, Float32.(e.val))
    end

    # Nucleus voxels (only if nuc_label is on).
    nuc_loc = Vector{Int32}[]
    nuc_val = Vector{Float32}[]
    if nuc_label
        for e in neur_vol.gp_nuc
            idxs, val = e
            push!(nuc_loc, Int32.(idxs))
            push!(nuc_val, fill(Float32(val), length(idxs)))
        end
    end

    freq_psf = psf_fft((H, W, D), psf; z_sub=scan_params.scan_avg)
    g = isnothing(g_blur) ? nothing : Matrix{Float32}(g_blur)
    return ScanVolume(soma_loc, soma_val, dend_loc, dend_val,
                      axon_loc, axon_val, nuc_loc, nuc_val,
                      freq_psf, (Np1, Np2, Np3), (H, W, D), g, nuc_label)
end

"""
    scan_volume_frame(scan_vol::ScanVolume, neur_act, scan_params::ScanParams;
                      z_off::Real=0, cutoff::Real=1e-2,
                      tpm_params::Union{Nothing,TPMParams}=nothing,
                      rng=Random.default_rng())
        -> (clean_img, z_off)

Synthesize one scan frame. `neur_act` is either a `NamedTuple` with
`soma`/`dend`/`bg` per-cell activity scalars (or vectors at the current
time step), or a single scalar applied uniformly. Returns the 2-D clean
image and the (possibly drifted) `z_off`. Ports `scan_volume_frame.m`.

If `tpm_params` is provided, the output is rescaled by the TPM
signal-scale via [`tpm_signal_scale`](@ref).
"""
function scan_volume_frame(scan_vol::ScanVolume,
                           neur_act,
                           scan_params::ScanParams;
                           z_off::Real=0,
                           cutoff::Real=1e-2,
                           tpm_params::Union{Nothing,TPMParams}=nothing,
                           rng::AbstractRNG=Random.default_rng())
    H, W, D = scan_vol.vol_sz
    Np1, Np2, Np3 = scan_vol.psf_sz
    sfrac = scan_params.sfrac

    # Optional small z-drift (motion model).
    z_off = scan_params.motion ?
        clamp(z_off + (rand(rng) < 0.005 ? 1 : 0) -
              (rand(rng) < 0.005 ? 1 : 0), -2, 2) :
        z_off
    z_loc = max(1, Int(floor(0.5 * (D - Np3))) + Int(round(z_off)))
    z_loc = min(z_loc, D - Np3 + 1)

    # Activity vectors.
    K_soma = length(scan_vol.soma_loc)
    K_axon = length(scan_vol.axon_loc)
    soma_act = _coerce_activity(neur_act, :soma, K_soma)
    dend_act = _coerce_activity(neur_act, :dend, K_soma)
    bg_act   = _coerce_activity(neur_act, :bg,   K_axon)
    nuc_act  = scan_vol.nuc_label ?
        _coerce_activity(neur_act, :soma, length(scan_vol.nuc_loc)) :
        Float32[]

    TMPvol = zeros(Float32, H, W, D)
    for ll in 1:K_soma
        a = soma_act[ll]
        a > cutoff || continue
        isempty(scan_vol.soma_loc[ll]) && continue
        for (i, li) in enumerate(scan_vol.soma_loc[ll])
            TMPvol[Int(li)] = scan_vol.soma_val[ll][i] * Float32(a)
        end
    end
    if scan_vol.nuc_label
        for ll in 1:length(scan_vol.nuc_loc)
            a = nuc_act[ll]
            a > cutoff || continue
            for (i, li) in enumerate(scan_vol.nuc_loc[ll])
                TMPvol[Int(li)] = scan_vol.nuc_val[ll][i] * Float32(a)
            end
        end
    end
    for ll in 1:K_soma
        a = dend_act[ll]
        a > cutoff || continue
        isempty(scan_vol.dend_loc[ll]) && continue
        for (i, li) in enumerate(scan_vol.dend_loc[ll])
            TMPvol[Int(li)] = scan_vol.dend_val[ll][i] * Float32(a)
        end
    end
    for ll in 1:K_axon
        a = bg_act[ll]
        a > cutoff || continue
        isempty(scan_vol.axon_loc[ll]) && continue
        for (i, li) in enumerate(scan_vol.axon_loc[ll])
            TMPvol[Int(li)] += scan_vol.axon_val[ll][i] * Float32(a)
        end
    end

    z_hi = min(z_loc + Np3 - 1, D)
    sub_vol = @view TMPvol[:, :, z_loc:z_hi]
    clean_img = (1 / (2 * sfrac^2)) .*
                single_scan(sub_vol, (Np1, Np2, Np3), scan_vol.freq_psf;
                            z_sub=scan_params.scan_avg, freq_opt=true)

    if !isnothing(scan_vol.g_blur)
        clean_img = clean_img .+ _conv2_same(clean_img, scan_vol.g_blur)
    end

    if sfrac > 1 && floor(sfrac) == sfrac
        # sum 2-D blocks of size (sfrac, sfrac), then subsample.
        sf = Int(sfrac)
        H2 = size(clean_img, 1) ÷ sf
        W2 = size(clean_img, 2) ÷ sf
        binned = zeros(Float32, H2, W2)
        for j in 1:W2, i in 1:H2
            s = 0f0
            for jj in 1:sf, ii in 1:sf
                s += clean_img[(i - 1) * sf + ii, (j - 1) * sf + jj]
            end
            binned[i, j] = s
        end
        clean_img = binned
    end

    if !isnothing(tpm_params)
        clean_img = clean_img .* Float32(tpm_signal_scale(tpm_params))
    end
    return clean_img, z_off
end

# Internal helper: pull a per-cell activity vector from `neur_act`,
# right-padded with zeros to length `K`. Vectors shorter than `K`
# correspond to "no activity for cells K_in+1 … K"; longer vectors are
# truncated.
function _coerce_activity(neur_act, field::Symbol, K::Integer)
    raw = if neur_act isa Number
        fill(Float32(neur_act), K)
    elseif neur_act isa AbstractVector
        Float32.(neur_act)
    elseif neur_act isa NamedTuple && hasproperty(neur_act, field)
        v = getproperty(neur_act, field)
        v isa Number ? fill(Float32(v), K) : Float32.(v)
    else
        zeros(Float32, K)
    end
    if length(raw) == K
        return raw
    elseif length(raw) < K
        padded = zeros(Float32, K)
        padded[1:length(raw)] .= raw
        return padded
    else
        return raw[1:K]
    end
end
