# Zernike polynomials and back-aperture aberration application.
#
# Derived from `zernike.m`, `generateZernike.m`, `applyZernike.m`, and
# `generateBA.m` in upstream NAOMi-Sim (Copyright 2021 Alex Song, Adam
# Charles, MIT). Zernike indexing follows the upstream `zidx` convention,
# which matches the standard Noll convention with `j = 1` → piston.

export zernike_polynomial, generate_zernike_weights, apply_zernike,
       generate_back_aperture

"""
    zernike_polynomial(j::Integer, x, y) -> Array

Evaluate the Noll-indexed Zernike polynomial `Z_j(x, y)` over a Cartesian
grid. `x` and `y` are coordinate arrays of identical shape, normalised
so the unit disk corresponds to `x² + y² ≤ 1`. The result has the same
shape as `x`.

Indexing follows upstream's `zidx`:

- `j = 1`: piston (Z = 1)
- `j = 2`: tilt-x (`2r cos θ`)
- `j = 3`: tilt-y (`2r sin θ`)
- `j = 4`: defocus (`√3 (2r² − 1)`)
- …

For all `j > 1`, the orthonormalisation constant is `√(2(n+1))` (or
`√(n+1)` for `m = 0`), so each `Z_j` integrates to `π` on the unit disk.
"""
function zernike_polynomial(j::Integer, x::AbstractArray, y::AbstractArray)
    size(x) == size(y) || throw(ArgumentError("x and y must have matching shapes"))
    n, m = _zidx(j)
    out = similar(x, Float64)
    for i in eachindex(x)
        xi, yi = x[i], y[i]
        r = sqrt(xi^2 + yi^2)
        if m == 0
            out[i] = sqrt(n + 1) * _zrf(n, 0, r)
        else
            θ = atan(yi, xi)
            base = sqrt(2 * (n + 1)) * _zrf(n, m, r)
            out[i] = iseven(j) ? base * cos(m * θ) : base * sin(m * θ)
        end
    end
    return out
end

# Upstream zidx -> (n, m). j is 1-based.
function _zidx(j::Integer)
    j ≥ 1 || throw(ArgumentError("Zernike index must be ≥ 1"))
    idx = ceil(Int, sqrt(0.25 + 2j) - 1.5)
    n = idx
    rem = j - (idx + 1) * idx ÷ 2
    m = if iseven(idx)
        2 * fld(rem, 2)
    else
        2 * cld(rem, 2) - 1
    end
    return n, m
end

# Radial polynomial R_n^m(r) (no orthonormalisation factor).
function _zrf(n::Integer, m::Integer, r::Real)
    R = 0.0
    half = (n - m) ÷ 2
    for s in 0:half
        sign = iseven(s) ? 1.0 : -1.0
        num = factorial(n - s)
        den = factorial(s) * factorial((n + m) ÷ 2 - s) * factorial((n - m) ÷ 2 - s)
        R += sign * (num / den) * r^(n - 2s)
    end
    return R
end

"""
    generate_zernike_weights(psf_params::PSFParams; offset = (0.0, 0.0))
        -> Vector{Float64}

Translate `psf_params.zernikeWt` into physical aberration amplitudes
(units of metres) by multiplying by `lambda * 1e-6` (upstream wavelength
is in microns, output is in metres). `offset` is reserved for the
`zernikeDst` spatial-modulation hook from upstream; not exercised here
(no `zernikeDst` field on `PSFParams`).
"""
function generate_zernike_weights(psf_params::PSFParams; offset = (0.0, 0.0))
    isempty(psf_params.zernikeWt) && return Float64[]
    return psf_params.zernikeWt .* (psf_params.lambda * 1e-6)
end

"""
    apply_zernike(Uin, X, Y, k, abb) -> Matrix{ComplexF64}

Apply a sum-of-Zernikes phase aberration to the scalar field `Uin`.
`X`, `Y` are normalised back-aperture coordinates (such that `X² + Y² ≤ 1`
inside the aperture); `k` is the wavenumber (rad/m); `abb` is a vector
of per-mode aberration amplitudes in metres (returned by
[`generate_zernike_weights`](@ref)).

The accumulated phase is mean-subtracted before being applied (matching
upstream — keeps the global phase offset out of the field).
"""
function apply_zernike(Uin::AbstractArray, X::AbstractArray, Y::AbstractArray,
                       k::Real, abb::AbstractVector)
    size(X) == size(Y) == size(Uin) ||
        throw(ArgumentError("Uin, X, Y must have matching shapes"))
    phase = zeros(Float64, size(X))
    for (j, a) in pairs(abb)
        a == 0 && continue
        phase .+= a .* zernike_polynomial(j, X, Y)
    end
    phase .-= sum(phase) / length(phase)
    return Uin .* exp.(1im .* k .* phase)
end

"""
    generate_back_aperture(vol_params::VolumeParams, psf_params::PSFParams;
                            vasc_sz = nothing) -> Matrix{ComplexF64}

Generate the Gaussian back-aperture field with the aberrations specified
by `psf_params.zernikeWt` applied. `vasc_sz` (the vasculature-bounding
box) is optional; if `nothing`, it is computed from
[`gaussian_beam_size`](@ref) the same way upstream `generateBA.m` does
when `vol_params.vasc_sz` is missing.

This is a port of the *simple* path of upstream `generateBA.m`
(`imax*jmax == 1`, no spatial Zernike distribution). The cell-array
path (one back-aperture per FOV pixel) is deferred.

Note upstream's `imax*jmax > 1` branch has a latent bug — it drops
the `X/objrad, Y/objrad` normalisation when calling `applyZernike`.
When that branch is ported we'll fix it.
"""
function generate_back_aperture(vol_params::VolumeParams, psf_params::PSFParams;
                                vasc_sz = nothing)
    if vasc_sz === nothing
        s = gaussian_beam_size(psf_params,
                               vol_params.vol_depth + vol_params.vol_sz[3] / 2)
        vasc_sz = (s[1] + vol_params.vol_sz[1],
                   s[2] + vol_params.vol_sz[2],
                   vol_params.vol_sz[3] + vol_params.vol_depth)
    end
    fl     = psf_params.obj_fl / 1000               # m
    D2     = 1e-6 * (1 / vol_params.vres) / psf_params.ss
    Nxy    = @. 1e-6 * (vasc_sz[1:2] - vol_params.vol_sz[1:2]) / D2
    Nxy    = (round(Int, Nxy[1]), round(Int, Nxy[2]))
    minN   = min(Nxy[1], Nxy[2])
    minN > 0 || throw(ArgumentError("Computed grid size is non-positive — check vol_params"))
    D1     = maximum(gaussian_beam_size(psf_params, fl * 1e6)) / minN / 1e6   # m
    nre    = psf_params.n
    rad    = tan(asin(psf_params.NA / nre)) * fl
    objrad = tan(asin(psf_params.objNA / nre)) * fl
    k      = 2 * nre * π / (psf_params.lambda * 1e-6)

    xs = ((-Nxy[1] ÷ 2):(Nxy[1] ÷ 2 - 1)) .* D1
    ys = ((-Nxy[2] ÷ 2):(Nxy[2] ÷ 2 - 1)) .* D1
    X = [x for x in xs, _ in ys]
    Y = [y for _ in xs, y in ys]
    Uout = generate_gaussian_profile(X, Y, rad, objrad, k, fl)
    abb = generate_zernike_weights(psf_params)
    return isempty(abb) ? Uout : apply_zernike(Uout, X ./ objrad, Y ./ objrad, k, abb)
end
