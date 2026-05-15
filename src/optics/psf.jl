# Gaussian-beam point-spread functions.
#
# Derived from `gaussian_psf.m`, `gaussian_psf_na.m`, `gaussianBeamSize.m`,
# and `generateGaussianProfile.m` in upstream NAOMi-Sim
# (Copyright 2021 Alex Song, Adam Charles, MIT). These functions are pure
# analytic evaluations on a grid; no FFTs (those appear in Chunk 13).

export gaussian_psf, gaussian_psf_na, gaussian_beam_size,
       generate_gaussian_profile

"""
    gaussian_psf(psflen, lambda, sampling, matSize; theta = 0)
        -> (psf, x, y, z)

Two-photon Gaussian point-spread function with an explicit axial half-power
length `psflen` (defined here as the plane where intensity drops to ½ of
the central plane, *not* the conventional FWHM). `lambda` is the
wavelength, `sampling` is `(sx, sy, sz)` (a 3-vector, or scalar broadcast
to all three), `matSize = (Nx, Ny, Nz)` is the output grid. `theta`
rotates the beam in the x–z plane (degrees).

Returns the squared intensity `psf` and the coordinate axes `x`, `y`, `z`
(centered such that `x[round(Nx/2)] = 0`).
"""
function gaussian_psf(psflen::Real, lambda::Real,
                      sampling::AbstractVector, matSize::AbstractVector;
                      theta::Real = 0)
    samp = length(sampling) < 3 ? fill(sampling[1], 3) : Float64.(sampling)
    Nx, Ny, Nz = matSize[1], matSize[2], matSize[3]
    n = 1.0
    zr = psflen / 2
    x = ((1:Nx) .- round(Int, Nx / 2)) .* samp[1]
    y = ((1:Ny) .- round(Int, Ny / 2)) .* samp[2]
    z = ((1:Nz) .- round(Int, Nz / 2)) .* samp[3]
    ct, st = cosd(theta), sind(theta)
    psf = Array{Float64}(undef, Nx, Ny, Nz)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        xv = x[i]; yv = y[j]; zv = z[k]
        xr =  ct * xv - st * zv
        zrot = st * xv + ct * zv
        denom = 1 + (zrot / zr)^2
        e = exp(-2 * π * n * (xr^2 + yv^2) / (zr * lambda * denom)) / denom
        psf[i, j, k] = e^2
    end
    return psf, collect(x), collect(y), collect(z)
end

gaussian_psf(psflen::Real, lambda::Real, sampling::Real,
             matSize::AbstractVector; kwargs...) =
    gaussian_psf(psflen, lambda, fill(sampling, 3), matSize; kwargs...)

"""
    gaussian_psf_na(na, lambda, sampling, matSize; theta = 0, nidx = 1.33)
        -> (psf, x, y, z, intensity)

NA-parameterised Gaussian PSF. `psflen` is computed as
`0.626 · lambda / (nidx − √(nidx² − na²))` (the standard 2-photon focal
extent); the rest of the formula matches [`gaussian_psf`](@ref) with `n`
set to `nidx`. Returns both the squared `psf` and the underlying complex
`intensity` field.

Note: this function's coordinate origin is at index `round(Nx/2) + 1`,
which is one sample off from [`gaussian_psf`](@ref) — a quirk of the
upstream definition (off-by-one in the `(0:matSize-1)` indexing); we
mirror it for fidelity.
"""
function gaussian_psf_na(na::Real, lambda::Real,
                         sampling::AbstractVector, matSize::AbstractVector;
                         theta::Real = 0, nidx::Real = 1.33)
    samp = length(sampling) < 3 ? fill(sampling[1], 3) : Float64.(sampling)
    Nx, Ny, Nz = matSize[1], matSize[2], matSize[3]
    psflen = 0.626 * lambda / (nidx - sqrt(nidx^2 - na^2))
    zr = psflen / 2
    x = ((0:Nx - 1) .- round(Int, Nx / 2)) .* samp[1]
    y = ((0:Ny - 1) .- round(Int, Ny / 2)) .* samp[2]
    z = ((0:Nz - 1) .- round(Int, Nz / 2)) .* samp[3]
    ct, st = cosd(theta), sind(theta)
    intensity = Array{Float64}(undef, Nx, Ny, Nz)
    psf = Array{Float64}(undef, Nx, Ny, Nz)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        xv = x[i]; yv = y[j]; zv = z[k]
        xr =  ct * xv - st * zv
        zrot = st * xv + ct * zv
        denom = 1 + (zrot / zr)^2
        e = exp(-2 * π * nidx * (xr^2 + yv^2) / (zr * lambda * denom)) / denom
        intensity[i, j, k] = e
        psf[i, j, k] = e^2
    end
    return psf, collect(x), collect(y), collect(z), intensity
end

gaussian_psf_na(na::Real, lambda::Real, sampling::Real,
                matSize::AbstractVector; kwargs...) =
    gaussian_psf_na(na, lambda, fill(sampling, 3), matSize; kwargs...)

"""
    gaussian_beam_size(psf_params::PSFParams, dist; apod = 2)
        -> (sx, sy, sz)

Conservative upper bound on the Gaussian beam waist at `dist` from focus,
returned as a 3-tuple `(sx, sx, 0)` (axial component is always zero — this
function is for lateral interaction-radius bookkeeping only).

`apod` is a linear scaling factor; the upstream default of 2 corresponds
to "about 5 apodizations" of margin.
"""
function gaussian_beam_size(psf_params::PSFParams, dist::Real; apod::Real = 2)
    s = ceil(Int, tan(asin(psf_params.objNA / psf_params.n)) * dist * 1.5) * apod
    return (s, s, 0)
end

"""
    generate_gaussian_profile(X, Y, rad, aper, k, fl, offset = (0, 0))
        -> Matrix{ComplexF64}

Gaussian back-aperture intensity profile with a fixed aperture and ideal
focusing phase. `X` and `Y` are coordinate matrices (broadcast-compatible
shapes; upstream uses 2-D `meshgrid` outputs). `rad` is the Gaussian 1/e²
radius, `aper` is the hard aperture radius, `k = 2π/λ` is the wavenumber,
`fl` is the focal length, and `offset = (dx, dy)` shifts the Gaussian
centre.
"""
function generate_gaussian_profile(X::AbstractArray, Y::AbstractArray,
                                   rad::Real, aper::Real, k::Real, fl::Real,
                                   offset = (0.0, 0.0))
    size(X) == size(Y) ||
        throw(ArgumentError("X and Y must have the same shape"))
    rho2 = X .^ 2 .+ Y .^ 2
    Uout = exp.(-((X .- offset[1]) .^ 2 .+ (Y .- offset[2]) .^ 2) ./ rad^2)
    Uout = Uout .* (rho2 .< aper^2)
    Uout = Uout .* exp.(-1im * k / (2 * fl) .* rho2)
    return Uout
end
