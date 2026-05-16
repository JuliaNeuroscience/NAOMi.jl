# Ideal-component spatial profiles + time-trace extraction.
#
# Ported from upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - calculateIdealComps.m  → calculate_ideal_comps
#   - comps2ideals.m         → comps2ideals
#   - times_from_profs.m     → times_from_profs
#
# `scan_ideal.m` references the unported (and absent in upstream)
# `single_scan_stack.m`; the working part of its functionality is
# subsumed by `calculate_ideal_comps`, so it is not separately ported.
# `constrainEstToSomas.m` post-processes a downstream-analysis `est`
# struct that this port has not constructed; deferred.

import Statistics
using LinearAlgebra: norm

export calculate_ideal_comps, comps2ideals, times_from_profs

"""
    calculate_ideal_comps(neur_vol, psf, neur_act, scan_params;
                          noise_params=nothing, tpm_params=nothing,
                          spike_opts=nothing, num_comps=nothing,
                          rng=Random.default_rng()) -> (comps, baseim, ideal)

Compute the per-cell spatial profile (`comps`), the baseline-activity
image (`baseim`), and the SNR-adjusted "ideal" profile stack (`ideal`).

`comps[:, :, k]` is the image that would be recorded if only cell `k`
were emitting fluorescence at its minimum-over-time level
`min(neur_act.soma[k, :])` (and likewise for its dendrite and the `k`-th
background process). Cells whose summed minimum activity is zero are
skipped — their slice in `comps` / `ideal` stays zero. Internally drives
[`scan_volume`](@ref) with `motion=false`.

Ports `calculateIdealComps.m`.
"""
function calculate_ideal_comps(neur_vol::NeuralVolume,
                               psf::AbstractArray{<:Real,3},
                               neur_act,
                               scan_params::ScanParams;
                               noise_params=nothing,
                               tpm_params=nothing,
                               spike_opts=nothing,
                               num_comps::Union{Nothing,Integer}=nothing,
                               rng::AbstractRNG=Random.default_rng())
    n0 = vec(minimum(neur_act.soma; dims=2))
    n1 = hasproperty(neur_act, :dend) ?
        vec(minimum(neur_act.dend; dims=2)) : copy(n0)
    n2 = hasproperty(neur_act, :bg) ?
        vec(minimum(neur_act.bg; dims=2)) : zeros(Float32, 1)
    K_soma = length(n0)
    K_axon = length(n2)

    g_idxs = falses(K_soma)
    K_min  = min(K_soma, length(n1), K_axon)
    for k in 1:K_min
        g_idxs[k] = (n0[k] + n1[k] + n2[k]) > 0
    end
    keep = findall(g_idxs)
    if !isnothing(num_comps)
        keep = keep[keep .<= num_comps]
    end
    n_keep = length(keep)
    n_keep == 0 && return (zeros(Float32, 0, 0, K_soma),
                           zeros(Float32, 0, 0),
                           zeros(Float32, 0, 0, K_soma))

    # Build diagonal "activate one component at its baseline level".
    soma_base = zeros(Float32, K_soma, n_keep)
    dend_base = zeros(Float32, K_soma, n_keep)
    bg_base   = zeros(Float32, K_axon, n_keep)
    for (j, k) in enumerate(keep)
        soma_base[k, j] = Float32(n0[k])
        dend_base[k, j] = Float32(n1[k])
        if k <= K_axon
            bg_base[k, j] = Float32(n2[k])
        end
    end

    sp_no_motion = ScanParams(; scan_buff=scan_params.scan_buff,
                              motion=false,
                              scan_avg=scan_params.scan_avg,
                              sfrac=scan_params.sfrac,
                              verbose=scan_params.verbose)

    base_act = (soma=soma_base, dend=dend_base, bg=bg_base)
    _, comps_kept = scan_volume(neur_vol, psf, base_act, sp_no_motion;
                                noise_params=noise_params,
                                tpm_params=tpm_params,
                                spike_opts=spike_opts,
                                rng=rng, return_clean=true)
    H, W = size(comps_kept, 1), size(comps_kept, 2)
    comps = zeros(Float32, H, W, K_soma)
    for (j, k) in enumerate(keep)
        comps[:, :, k] .= @view comps_kept[:, :, j]
    end

    # Baseline image: scan with everyone at their minimum simultaneously.
    f0_act = (soma=reshape(Float32.(n0), :, 1),
              dend=reshape(Float32.(n1), :, 1),
              bg=reshape(Float32.(n2), :, 1))
    _, base_mov = scan_volume(neur_vol, psf, f0_act, sp_no_motion;
                              noise_params=noise_params,
                              tpm_params=tpm_params,
                              spike_opts=spike_opts,
                              rng=rng, return_clean=true)
    baseim = @view base_mov[:, :, 1]

    ideal = comps2ideals(comps, Array(baseim))
    return comps, Array(baseim), ideal
end

"""
    comps2ideals(comps::AbstractArray{<:Real,3},
                 baseim::AbstractMatrix{<:Real}; k::Real=2) -> Array

SNR-threshold each `comps[:, :, i]` against `baseim`. Per component, the
ratio image `r = comps[:, :, i] ./ baseim` is computed, a cutoff is
derived from the mean of the top-5 values of `r`, the largest connected
component above that cutoff is kept (if at least 5 pixels), and the
masked `comps[:, :, i]` is returned. Ports `comps2ideals.m`.
"""
function comps2ideals(comps::AbstractArray{<:Real,3},
                      baseim::AbstractMatrix{<:Real}; k::Real=2)
    H, W, K = size(comps)
    min_num_el = 5
    ideal = similar(Array(comps), Float32)
    fill!(ideal, 0)
    # Ratio image, with division-by-zero giving Inf/NaN we treat as zero.
    for c in 1:K
        ratio = Vector{Float32}(undef, H * W)
        for j in 1:W, i in 1:H
            b = Float32(baseim[i, j])
            v = b == 0 ? 0f0 : Float32(comps[i, j, c]) / b
            ratio[(j - 1) * H + i] = isfinite(v) ? v : 0f0
        end
        # cutoff from top-5 mean.
        finite_ratio = filter(isfinite, ratio)
        isempty(finite_ratio) && continue
        sorted = sort(finite_ratio; rev=true)
        top_mean = Statistics.mean(sorted[1:min(min_num_el, length(sorted))])
        top_mean > 0 || continue
        cutoff = 1 / (k + 2 / top_mean)
        mask = reshape(ratio .> cutoff, H, W)
        any(mask) || continue
        comp_id, sizes = _label_connected(mask)
        isempty(sizes) && continue
        biggest = argmax(sizes)
        sizes[biggest] >= min_num_el || continue
        for j in 1:W, i in 1:H
            if comp_id[i, j] == biggest
                ideal[i, j, c] = Float32(comps[i, j, c])
            end
        end
    end
    return ideal
end

# 4-connectivity connected-component labelling on a binary mask. Returns
# the labelled array (0 = background) and the per-label pixel count.
function _label_connected(mask::AbstractMatrix{Bool})
    H, W = size(mask)
    labels = zeros(Int, H, W)
    sizes = Int[]
    next = 0
    stack = Tuple{Int,Int}[]
    for j in 1:W, i in 1:H
        (mask[i, j] && labels[i, j] == 0) || continue
        next += 1
        push!(sizes, 0)
        push!(stack, (i, j))
        labels[i, j] = next
        while !isempty(stack)
            (ii, jj) = pop!(stack)
            sizes[next] += 1
            for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))
                ni, nj = ii + di, jj + dj
                (1 <= ni <= H && 1 <= nj <= W) || continue
                if mask[ni, nj] && labels[ni, nj] == 0
                    labels[ni, nj] = next
                    push!(stack, (ni, nj))
                end
            end
        end
    end
    return labels, sizes
end

"""
    times_from_profs(mov::AbstractArray{<:Real,3},
                     neur_prof::AbstractArray{<:Real,3};
                     bg_profs=nothing, lambda::Real=0,
                     nnls::Bool=true, max_iter::Integer=200,
                     tol::Real=1e-8)
        -> (x_est, x_bg)

Recover per-component activity traces from a movie `mov` (H × W × T) given
a stack of spatial profiles `neur_prof` (H × W × C). `bg_profs` is an
optional H × W × B background profile stack; if `nothing`, the median
frame of `mov` is used (normalised to unit Frobenius norm). With
`lambda == 0` and `nnls=true` the per-time-step problem
`argmin_{x ≥ 0} ‖A x − y‖²` is solved by projected gradient; with
`nnls=false` it falls back to unconstrained least squares. The L1-
penalised (`lambda > 0`) and TFOCS-specific paths from upstream are
deferred.

Returns `(x_est, x_bg)` where `x_est` is `C × T` and `x_bg` is `B × T`.
Ports the default branches of `times_from_profs.m`.
"""
function times_from_profs(mov::AbstractArray{<:Real,3},
                          neur_prof::AbstractArray{<:Real,3};
                          bg_profs::Union{Nothing,AbstractArray{<:Real,3}}=nothing,
                          lambda::Real=0,
                          nnls::Bool=true,
                          max_iter::Integer=200,
                          tol::Real=1e-8)
    lambda > 0 && error("L1-penalised path (lambda > 0) not ported; pass lambda=0")
    H, W, T = size(mov)
    C = size(neur_prof, 3)
    (size(neur_prof, 1) == H && size(neur_prof, 2) == W) ||
        throw(DimensionMismatch("neur_prof and mov must share the lateral extent"))

    bg = if isnothing(bg_profs)
        med = Array{Float32}(undef, H, W)
        # Per-pixel median over time.
        buf = Vector{Float32}(undef, T)
        for j in 1:W, i in 1:H
            for t in 1:T
                buf[t] = Float32(mov[i, j, t])
            end
            med[i, j] = Statistics.median(buf)
        end
        nm = sqrt(sum(abs2, med))
        nm > 0 ? reshape(med ./ nm, H, W, 1) : reshape(med, H, W, 1)
    else
        Float32.(bg_profs)
    end
    B = size(bg, 3)

    A = zeros(Float32, H * W, C + B)
    for c in 1:C
        for j in 1:W, i in 1:H
            A[(j - 1) * H + i, c] = Float32(neur_prof[i, j, c])
        end
    end
    for b in 1:B
        for j in 1:W, i in 1:H
            A[(j - 1) * H + i, C + b] = Float32(bg[i, j, b])
        end
    end

    Y = Matrix{Float32}(undef, H * W, T)
    for t in 1:T, j in 1:W, i in 1:H
        Y[(j - 1) * H + i, t] = Float32(mov[i, j, t])
    end

    X = if !nnls
        max.(Float32.(A \ Y), 0f0)
    else
        _projected_gradient_nnls(A, Y; max_iter=max_iter, tol=tol)
    end

    x_est = X[1:C, :]
    x_bg  = X[(C + 1):(C + B), :]
    return x_est, x_bg
end

# Projected gradient NNLS, columns of Y solved jointly.
# x_{k+1} = max(0, x_k - α (AᵀA x_k - AᵀY)), α = 1 / ‖AᵀA‖.
function _projected_gradient_nnls(A::AbstractMatrix{<:Real},
                                  Y::AbstractMatrix{<:Real};
                                  max_iter::Integer=200, tol::Real=1e-8)
    AtA = Float32.(A' * A)
    AtY = Float32.(A' * Y)
    # Spectral-norm upper bound via power iteration.
    n = size(AtA, 1)
    v = randn(Float32, n)
    v ./= max(norm(v), eps(Float32))
    local lam
    for _ in 1:30
        w = AtA * v
        nw = norm(w)
        nw == 0 && (lam = 1f0; break)
        v = w ./ nw
        lam = nw
    end
    lam = max(lam, eps(Float32))
    α = 1 / Float32(lam)
    X = max.(AtY ./ Float32(lam), 0f0)
    prev = copy(X)
    for it in 1:max_iter
        grad = AtA * X .- AtY
        X .= max.(X .- α .* grad, 0f0)
        if it % 5 == 0
            δ = maximum(abs.(X .- prev)) / (maximum(abs.(X)) + eps(Float32))
            δ < tol && return X
            prev .= X
        end
    end
    return X
end
