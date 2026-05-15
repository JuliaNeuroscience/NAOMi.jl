# Fresnel propagation, z-projection, beam-width measurement, two-photon
# photon-flux estimation, and collection-side hemoglobin absorption.
#
# Derived from `fresnel_propagation_multi.m`, `groupzproject.m`,
# `widthestimate.m`, `widthestimate3D.m`, `tpmSignalscale.m`, and the
# collection-mask block of `simulate_optical_propagation.m` in upstream
# NAOMi-Sim (Copyright 2021 Alex Song, Adam Charles, MIT).
#
# The full per-tile cortical-light-path orchestrators
# (`genCorticalLightPath.m`, `genCorticalLightPathLite.m`) and the
# top-level `simulate_optical_propagation.m` driver are deferred to a
# later chunk; see ANALYSIS_PLAN.md.

export fresnel_propagation_multi, group_z_project,
       width_estimate, width_estimate_3d,
       tpm_signal_scale, collection_mask

"""
    fresnel_propagation_multi(Uin, lambda, dx, z, phi, nidx; saveall = false)
        -> Uout  (or (Uout, UoutAll) when `saveall = true`)

Multi-step angular-spectrum propagation of a complex scalar field through a
sequence of phase-modulated planes. Inputs:

- `Uin`     — incident field `(Nx, Ny)` complex matrix.
- `lambda`  — vacuum wavelength (m).
- `dx`      — length-`n` vector of lateral pixel pitches (m), one per step.
- `z`       — length-`n` vector of axial positions (m); `z[1]` is the
              source plane and `z[end]` the final plane.
- `phi`     — `(Nx, Ny, m)` complex phase modulations applied at each
              plane. When `m == n` the modulation at each plane (including
              both endpoints) is applied; when `m == n − 1` the final plane
              skips modulation (matches upstream's `size(phi,3) < n`
              special case).
- `nidx`    — refractive index of the propagation medium.

The kernel matches `fresnel_propagation_multi.m` from upstream NAOMi-Sim:
quadratic phase ramp at entry, frequency-domain quadratic transfer
function `Q` per step, and a final quadratic correction at exit.

Returns the field after the final plane. When `saveall = true`, also
returns the per-plane stack `UoutAll`.
"""
function fresnel_propagation_multi(Uin::AbstractMatrix{<:Number},
                                   lambda::Real,
                                   dx::AbstractVector,
                                   z::AbstractVector,
                                   phi::AbstractArray{<:Number,3},
                                   nidx::Real;
                                   saveall::Bool = false)
    n = length(z)
    n >= 2 || throw(ArgumentError("z must have at least two elements"))
    length(dx) == n ||
        throw(ArgumentError("dx must have the same length as z"))
    Nx, Ny = size(Uin)
    (size(phi, 1), size(phi, 2)) == (Nx, Ny) ||
        throw(ArgumentError("phi must have the same lateral shape as Uin"))
    mphi = size(phi, 3)
    (mphi == n || mphi == n - 1) ||
        throw(ArgumentError("size(phi,3) must equal length(z) or length(z)-1"))

    λ = lambda / nidx
    k = 2π / λ
    df = 1.0 ./ (Nx .* dx)
    dz = diff(z)
    sc = dx[2:end] ./ dx[1:end-1]

    nxs = collect(-Nx ÷ 2 : Nx ÷ 2 - 1)
    nys = collect(-Ny ÷ 2 : Ny ÷ 2 - 1)
    nx = [Float64(x) for x in nxs, _ in nys]
    ny = [Float64(y) for _ in nxs, y in nys]

    # Apply entry quadratic phase + first plane modulation.
    U = ComplexF64.(Uin) .*
        exp.(1im .* k .* ((nx .* dx[1]) .^ 2 .+ (ny .* dx[1]) .^ 2) .*
             (1 - sc[1]) ./ (2 * dz[1])) .*
        @view(phi[:, :, 1])

    Uall = saveall ? Array{ComplexF64}(undef, Nx, Ny, n) : nothing
    if saveall
        @views Uall[:, :, 1] .= U
    end

    fft_plan  = plan_fft(U)
    ifft_plan = plan_ifft(U)

    tol = 1e-12
    Q = exp.(-1im .* π .* λ .* dz[1] .*
             ((nx .* df[1]) .^ 2 .+ (ny .* df[1]) .^ 2) ./ sc[1])

    for i in 1:(n - 1)
        if i > 1 &&
           !(abs(dz[i] - dz[i - 1]) < tol &&
             abs(df[i] - df[i - 1]) < tol &&
             abs(sc[i] - sc[i - 1]) < tol)
            Q = exp.(-1im .* π .* λ .* dz[i] .*
                     ((nx .* df[i]) .^ 2 .+ (ny .* df[i]) .^ 2) ./ sc[i])
        end

        Uf = fftshift(fft_plan * fftshift(U ./ sc[i]))
        Uf .= Q .* Uf
        Up = ifftshift(ifft_plan * ifftshift(Uf))

        if i == n - 1 && mphi < n
            U = Up
        else
            U = @view(phi[:, :, i + 1]) .* Up
        end

        if saveall
            @views Uall[:, :, i + 1] .= U
        end
    end

    Uout = U .* exp.(1im .* k ./ 2 .* (sc[end] - 1) ./ (sc[end] .* dz[end]) .*
                     ((nx .* dx[end]) .^ 2 .+ (ny .* dx[end]) .^ 2))

    return saveall ? (Uout, Uall) : Uout
end

# Convenience: scalar `dx` broadcast to all `n` planes.
function fresnel_propagation_multi(Uin::AbstractMatrix{<:Number},
                                   lambda::Real, dx::Real,
                                   z::AbstractVector,
                                   phi::AbstractArray{<:Number,3},
                                   nidx::Real; kwargs...)
    return fresnel_propagation_multi(Uin, lambda, fill(Float64(dx), length(z)),
                                     z, phi, nidx; kwargs...)
end

"""
    group_z_project(images, groupsize; type = :mean) -> Array

Project a 3-D array along its third dimension in groups of `groupsize`.
`type` is one of `:sum`, `:prod`, `:mean`, `:max`, `:min`. If
`size(images, 3)` is not divisible by `groupsize`, the final group
contains the remainder (matching upstream `groupzproject.m`).
"""
function group_z_project(images::AbstractArray{<:Any,3}, groupsize::Integer;
                         type::Symbol = :mean)
    op = if type === :sum
        sum
    elseif type === :prod
        prod
    elseif type === :mean
        a -> sum(a) / length(a)
    elseif type === :max
        maximum
    elseif type === :min
        minimum
    else
        throw(ArgumentError("Unsupported projection type :$type " *
                            "(supported: :sum, :prod, :mean, :max, :min)"))
    end

    Nx, Ny, Nz = size(images)
    ng_full = Nz ÷ groupsize
    rest = Nz - ng_full * groupsize
    ng = rest == 0 ? ng_full : ng_full + 1

    T = eltype(images)
    Tout = T <: Integer ? Float64 : float(T)
    out = Array{Tout}(undef, Nx, Ny, ng)

    @inbounds for g in 1:ng_full
        lo = (g - 1) * groupsize + 1
        hi = g * groupsize
        for j in 1:Ny, i in 1:Nx
            out[i, j, g] = op(@view(images[i, j, lo:hi]))
        end
    end
    if rest > 0
        lo = Nz - rest + 1
        @inbounds for j in 1:Ny, i in 1:Nx
            out[i, j, end] = op(@view(images[i, j, lo:Nz]))
        end
    end
    return out
end

"""
    width_estimate(vector, fraction = 0.5) -> Real

Continuous-domain width of a 1-D peak at the level `fraction` of its
maximum, found by linear interpolation between the two threshold
crossings. With `fraction = 0.5` this is the FWHM. The signal must
cross the threshold exactly twice (a simple unimodal peak); upstream
`widthestimate.m` makes the same assumption.
"""
function width_estimate(vector::AbstractVector{<:Real}, fraction::Real = 0.5)
    v = vector ./ maximum(vector)
    greater = v .> fraction
    flips = findall(d -> d != 0, diff(Int.(greater)))
    length(flips) >= 2 ||
        throw(ArgumentError("could not find two threshold crossings " *
                            "(signal must cross fraction exactly twice)"))
    f1, f2 = flips[1], flips[2]
    s1 = (v[f1 + 1] - fraction) / (v[f1 + 1] - v[f1])
    s2 = (v[f2]     - fraction) / (v[f2]     - v[f2 + 1])
    # Width between linear-interp threshold crossings at positions
    # `(f1 + 1 − s1)` and `(f2 + s2)` is `(f2 − f1 − 1) + s1 + s2`.
    # Upstream `widthestimate.m` returns `(f2 − f1) + s1 + s2`, which is
    # one sample too wide; we correct it here.
    return s1 + s2 + (f2 - f1 - 1)
end

"""
    width_estimate_3d(matrix, fraction = 0.5) -> NTuple{3,Real}

Apply [`width_estimate`](@ref) to each of the three lateral/axial
sum-projections of a 3-D peak. Returns `(width_x, width_y, width_z)`.
"""
function width_estimate_3d(matrix::AbstractArray{<:Real,3},
                           fraction::Real = 0.5)
    v1 = vec(sum(matrix, dims = (2, 3)))
    v2 = vec(sum(matrix, dims = (1, 3)))
    v3 = vec(sum(matrix, dims = (1, 2)))
    return (width_estimate(v1, fraction),
            width_estimate(v2, fraction),
            width_estimate(v3, fraction))
end

"""
    tpm_signal_scale(tpm_params::TPMParams,
                      psf_params::Union{PSFParams,Nothing} = nothing)
        -> Float64

Average number of fluorescence photons collected per second under the
two-photon imaging model of Xu and Webb (1996, JOSA B). Equation:

```
<F(t)> = (1/2) · phi · eta · C · delta · gp · 8 · nidx · <P(t)>² /
         (f · tau · pi · lambda)
```

When `psf_params` is supplied, `nidx`, `nac`, and `lambda` are taken
from it rather than from `tpm_params` (matching upstream
`tpmSignalscale.m`'s two-argument signature).
"""
function tpm_signal_scale(tpm_params::TPMParams,
                          psf_params::Union{PSFParams,Nothing} = nothing)
    if psf_params === nothing
        nidx   = tpm_params.nidx
        lambda = tpm_params.lambda
    else
        nidx   = psf_params.n
        lambda = psf_params.lambda
    end

    phi   = tpm_params.phi
    eta   = tpm_params.eta
    conc  = tpm_params.conc * 1e-6 * 6.02e23 * 1e3
    delta = tpm_params.delta * 1e-58
    gp    = tpm_params.gp
    f     = tpm_params.f * 1e6
    tau   = tpm_params.tau * 1e-15
    pavg  = tpm_params.pavg
    λm    = lambda * 1e-6
    # Power → photons per second: P / (h c / λ) with P in W.
    pphot = 1e-3 * pavg / (6.626e-34 * 3e8 / λm)
    return phi * eta * conc * delta * gp * 8 * nidx * pphot^2 /
           (2 * f * tau * π * λm)
end

"""
    collection_mask(vol_params::VolumeParams, psf_params::PSFParams,
                     neur_ves; vasc_sz = nothing) -> Matrix{Float64}

2-D map of two-photon emission collection efficiency through the
vasculature, modelled by Beer–Lambert absorption integrated over the
collection cone (Equation 4-equivalent of upstream
`simulate_optical_propagation.m`, lines 415–442).

Inputs:

- `neur_ves` — 3-D vasculature volume (`vasc_sz · vres` in shape) with
  voxel values in `[0, 1]` denoting fractional blood occupancy.
- `vasc_sz`  — optional explicit `(sx, sy, sz)` bounding box in µm; if
  omitted it is derived from [`gaussian_beam_size`](@ref).

Returns a `(vol_sz[1]·vres, vol_sz[2]·vres)` matrix of transmission
factors `T ∈ [0, 1]` (`T = 10^(−ℓ · hemoabs / vres)`, where `ℓ` is the
column-integrated vessel path within the collection cone).

Notes:

- The structuring element at each depth is a disk of radius
  `tan(asin(objNA/n)) · depth · vres`, normalised to unit sum, then
  convolved (with zero-padded boundaries) against the `proppx`-block
  sum-projection of `neur_ves`.
- The aperture-cropping autocrop step in upstream is a no-op for the
  result; it only trims the kernel for efficiency.
"""
function collection_mask(vol_params::VolumeParams, psf_params::PSFParams,
                         neur_ves::AbstractArray{<:Real,3};
                         vasc_sz = nothing)
    vres      = vol_params.vres
    vol_sz    = vol_params.vol_sz
    vol_depth = vol_params.vol_depth

    if vasc_sz === nothing
        s = gaussian_beam_size(psf_params, vol_depth + vol_sz[3] / 2)
        vasc_sz = (s[1] + vol_sz[1], s[2] + vol_sz[2],
                   vol_sz[3] + vol_depth)
    end

    proppx_f = psf_params.prop_sz * vres
    proppx   = round(Int, proppx_f)
    proppx > 0 || throw(ArgumentError("prop_sz × vres must be ≥ 1"))

    vascpx = (round(Int, vasc_sz[1] * vres),
              round(Int, vasc_sz[2] * vres),
              round(Int, vasc_sz[3] * vres))
    volpx  = (round(Int, vol_sz[1] * vres),
              round(Int, vol_sz[2] * vres),
              round(Int, vol_sz[3] * vres))

    # Centre the (smaller or equal) input vasculature in the vasc-sized array
    # before z-projecting. Upstream's full-pipeline path arranges this in
    # `simulate_optical_propagation.m`; here we accept either an
    # already-vasc-sized or vol-sized neur_ves.
    nv = if size(neur_ves) == vascpx
        neur_ves
    elseif size(neur_ves)[1:2] == volpx[1:2]
        padded = zeros(eltype(neur_ves), vascpx[1], vascpx[2], size(neur_ves, 3))
        dx_off = (vascpx[1] - volpx[1]) ÷ 2
        dy_off = (vascpx[2] - volpx[2]) ÷ 2
        padded[dx_off .+ (1:volpx[1]), dy_off .+ (1:volpx[2]), :] .= neur_ves
        padded
    else
        throw(ArgumentError("neur_ves shape $(size(neur_ves)) is neither " *
                            "vasc-sized $(vascpx) nor vol-sized $(volpx)"))
    end

    tmp = group_z_project(Float64.(nv), proppx; type = :sum)
    nz_blocks = size(tmp, 3)

    coldist = vres .* tan(asin(psf_params.objNA / psf_params.n)) .*
              ((vol_depth + vol_sz[3] / 2) .-
               (proppx / vres .* ((1:nz_blocks) .- 0.5)))

    colmask = zeros(Float64, vascpx[1], vascpx[2])
    for iz in 1:nz_blocks
        coldist[iz] > 0 || continue
        K = _disk_kernel(coldist[iz])
        colmask .+= _conv2_same(@view(tmp[:, :, iz]), K)
    end

    # Crop to volume footprint.
    dx_off = (vascpx[1] - volpx[1]) ÷ 2
    dy_off = (vascpx[2] - volpx[2]) ÷ 2
    cropped = colmask[dx_off .+ (1:volpx[1]), dy_off .+ (1:volpx[2])]
    return 10.0 .^ (-cropped ./ vres .* psf_params.hemoabs)
end

# Disk kernel of given radius (in pixels), normalised to unit sum.
function _disk_kernel(radius::Real)
    r = max(0, ceil(Int, radius))
    sz = 2r + 1
    K = zeros(Float64, sz, sz)
    @inbounds for j in 1:sz, i in 1:sz
        if (i - r - 1)^2 + (j - r - 1)^2 <= radius^2
            K[i, j] = 1.0
        end
    end
    s = sum(K)
    return s > 0 ? K ./ s : K
end

# 2-D zero-padded correlation matching MATLAB's `conv2(A, K, 'same')`
# for a symmetric kernel (centre-aligned, zero outside-array support).
function _conv2_same(A::AbstractMatrix{<:Real}, K::AbstractMatrix{<:Real})
    Nx, Ny = size(A)
    kx, ky = size(K)
    cx, cy = (kx + 1) ÷ 2, (ky + 1) ÷ 2
    out = zeros(Float64, Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        s = 0.0
        for jk in 1:ky, ik in 1:kx
            ia = i + ik - cx
            ja = j + jk - cy
            (1 <= ia <= Nx && 1 <= ja <= Ny) || continue
            s += A[ia, ja] * K[ik, jk]
        end
        out[i, j] = s
    end
    return out
end
