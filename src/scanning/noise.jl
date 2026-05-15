# Measurement-noise model.
#
# Ported from upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - PoissonGaussNoiseModel.m → poisson_gauss_noise
#   - applyNoiseModel.m        → apply_noise_model
#   - pixel_bleed.m            → pixel_bleed

export poisson_gauss_noise, pixel_bleed, apply_noise_model

"""
    poisson_gauss_noise(clean_in::AbstractMatrix{<:Real},
                        noise_params::NoiseParams;
                        rng=Random.default_rng()) -> Matrix{Float32}

Apply the Poisson-then-lognormal-then-Gaussian noise model:

```
count   ~ Poisson(clean_in + darkcount)
m       = count · μ
v       = count · σ
log_μ   = log(m² / √(v + m²))
log_σ   = √(log(v / m² + 1))
intermediate ~ LogNormal(log_μ, log_σ)
noisy   = round(intermediate + N(μ₀, σ₀²))
```

NaNs (from `count == 0`) are mapped to `0`. Ports
`PoissonGaussNoiseModel.m`.
"""
function poisson_gauss_noise(clean_in::AbstractMatrix{<:Real},
                             noise_params::NoiseParams;
                             rng::AbstractRNG=Random.default_rng())
    H, W = size(clean_in)
    out  = Matrix{Float32}(undef, H, W)
    @inbounds for j in 1:W, i in 1:H
        lam = max(0.0f0, Float32(clean_in[i, j]) + Float32(noise_params.darkcount))
        cnt = Float32(rand(rng, Poisson(lam)))
        m = cnt * Float32(noise_params.mu)
        v = cnt * Float32(noise_params.sigma)
        if cnt == 0
            ln_val = 0f0
        else
            μ2 = log(m^2 / sqrt(v + m^2))
            σ2 = sqrt(log(v / m^2 + 1))
            ln_val = Float32(exp(μ2 + σ2 * randn(rng)))
            isnan(ln_val) && (ln_val = 0f0)
        end
        gauss = Float32(noise_params.mu0) + Float32(noise_params.sigma0) * randn(rng, Float32)
        out[i, j] = round(ln_val + gauss)
    end
    return out
end

"""
    pixel_bleed(frame::AbstractMatrix{<:Real}, p::Real, b_max::Real;
                rng=Random.default_rng()) -> Matrix{Float32}

Model per-pixel electronic bleed-through: with probability `p` each
pixel donates a fraction (uniform in `(0, b_max]`) of its value to the
*next* pixel in scan order (row-major scan: same row, next column;
end-of-row wraps to next row's first column). Ports `pixel_bleed.m`.
"""
function pixel_bleed(frame::AbstractMatrix{<:Real}, p::Real, b_max::Real;
                     rng::AbstractRNG=Random.default_rng())
    H, W = size(frame)
    fr = Float32.(frame)
    if p <= 0
        return fr
    end
    # x_bleed = b_max * max(rand - (1-p), 0) / p — Bernoulli-gated uniform.
    x_bleed = Matrix{Float32}(undef, H, W)
    @inbounds for j in 1:W, i in 1:H
        r = rand(rng, Float32)
        x_bleed[i, j] = max(r - Float32(1 - p), 0f0) * Float32(b_max / p)
    end
    # Shifted versions (the value at (i,j) bleeds INTO (i,j+1) or
    # (i+1, 1) if at row end). The MATLAB expression
    # `[[0; x_bleed(1:end-1, end)], x_bleed(:, 1:end-1)]` constructs the
    # previous-pixel matrix: column 1 is the previous frame's last column
    # shifted down, columns 2..end are x_bleed columns 1..end-1.
    prev_xb = Matrix{Float32}(undef, H, W)
    prev_fr = Matrix{Float32}(undef, H, W)
    @inbounds prev_xb[1, 1] = 0f0
    @inbounds prev_fr[1, 1] = 0f0
    @inbounds for i in 2:H
        prev_xb[i, 1] = x_bleed[i - 1, end]
        prev_fr[i, 1] = fr[i - 1, end]
    end
    @inbounds for j in 2:W, i in 1:H
        prev_xb[i, j] = x_bleed[i, j - 1]
        prev_fr[i, j] = fr[i, j - 1]
    end
    out = fr .- x_bleed .* fr .+ prev_xb .* prev_fr
    return out
end

"""
    apply_noise_model(clean_mov::AbstractArray{<:Real,3}, noise_params::NoiseParams;
                      rng=Random.default_rng()) -> Array{Float32,3}

Apply [`poisson_gauss_noise`](@ref) frame-by-frame to a movie, followed
by [`pixel_bleed`](@ref). Returns a `Float32` movie of the same size.
Ports `applyNoiseModel.m` (Poisson-Gauss branch only; the dynode-chain
branch is upstream-specific and not ported).
"""
function apply_noise_model(clean_mov::AbstractArray{<:Real,3},
                           noise_params::NoiseParams;
                           rng::AbstractRNG=Random.default_rng())
    H, W, T = size(clean_mov)
    mov = Array{Float32, 3}(undef, H, W, T)
    @inbounds for k in 1:T
        noisy = poisson_gauss_noise(@view(clean_mov[:, :, k]), noise_params; rng=rng)
        mov[:, :, k] .= pixel_bleed(noisy, noise_params.bleedp,
                                    noise_params.bleedw; rng=rng)
    end
    return mov
end
