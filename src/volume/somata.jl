# Soma generation.
#
# Ported from the upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - generateNeuralBody.m         → generate_neural_body
#   - sampleDenseNeurons.m         → sample_dense_neurons
#   - generateNeuralVolume.m       → generate_neural_volume
#   - isolateVisibleSomas.m        → isolate_visible_somas
#   - teardrop_poj.m               → teardrop_projection
#   - SpiralSampleSphere.m         → spiral_sample_sphere
#   - masked_3DGP_v2.m             → masked_3d_gp

using LinearAlgebra

export spiral_sample_sphere, teardrop_projection, generate_neural_body,
       sample_dense_neurons, generate_neural_volume,
       isolate_visible_somas, point_in_soma, masked_3d_gp

# ---------------------------------------------------------------------------
# Sphere sampling and triangulation-free shape primitives
# ---------------------------------------------------------------------------

"""
    spiral_sample_sphere(N::Integer) -> Matrix{Float64}

Generate `N` approximately-uniform samples on the unit sphere using the
spiral / golden-angle method (Carlson, 2011). Returns an `N × 3` matrix
whose rows are unit vectors.

Ports the upstream `SpiralSampleSphere.m` from the
`S2_Sampling_Suite` bundle, without triangulation — downstream code in
this port uses a star-shape radial test (`point_in_soma`) rather than
mesh-based `intriangulation`.
"""
function spiral_sample_sphere(N::Integer)
    N ≥ 2 || throw(ArgumentError("N must be ≥ 2"))
    gr = (1 + sqrt(5)) / 2
    ga = 2π * (1 - 1 / gr)
    V = Matrix{Float64}(undef, N, 3)
    for i in 0:(N - 1)
        lat = acos(1 - 2 * i / (N - 1))
        lon = i * ga
        V[i + 1, 1] = sin(lat) * cos(lon)
        V[i + 1, 2] = sin(lat) * sin(lon)
        V[i + 1, 3] = cos(lat)
    end
    return V
end

"""
    teardrop_projection(V::AbstractMatrix, p::Real=1) -> Matrix

Project unit-sphere samples to a teardrop (pyramidal-cell mean shape),
following upstream `teardrop_poj.m`. For each input row `(x, y, z)`,

```
r  = sqrt(x² + y²)
θ  = π − atan2(r, |z|) − (z>0)·π
out = ( (x/r)·sin(θ)·sin(θ/2)^p,
        (y/r)·sin(θ)·sin(θ/2)^p,
        −cos(θ) )
```

NaNs from `0/0` near the poles are mapped to `0`.
"""
function teardrop_projection(V::AbstractMatrix{<:Real}, p::Real=1)
    N = size(V, 1)
    out = zeros(Float64, N, 3)
    for i in 1:N
        x, y, z = V[i, 1], V[i, 2], V[i, 3]
        rr = sqrt(x^2 + y^2)
        tt = π - atan(rr / abs(z)) - (z > 0 ? π : 0.0)
        # the 0.5 coefficient inside sin matches upstream's commented-active line.
        if rr > 0
            s = sin(tt) * sin(0.5 * tt)^p
            out[i, 1] = (x / rr) * s
            out[i, 2] = (y / rr) * s
        end
        out[i, 3] = -cos(tt)
    end
    return out
end

# ---------------------------------------------------------------------------
# Generate a single neural body (soma + nucleus mesh)
# ---------------------------------------------------------------------------

"""
    generate_neural_body(neur_params::NeuronParams;
                         V=nothing, dists=nothing, Rtear=nothing,
                         rng=Random.default_rng())
        -> (Vcell, Vnuc, rotation)

Sample a single neuron's surface mesh from an isotropic Gaussian process
on the sphere, optionally biased by a teardrop mean (for
`neur_type === :pyr`). Returns:

- `Vcell` — `N × 3` matrix of soma surface vertices.
- `Vnuc`  — `N × 3` matrix of nucleus surface vertices.
- `rotation` — 3-element vector of Euler rotation angles (degrees) used
  to rotate both meshes.

The optional `V` argument lets callers pass a precomputed spiral sphere
sampling (reused across many cells). `dists` is an optional precomputed
matrix of pairwise sphere distances. Ports `generateNeuralBody.m`.
"""
function generate_neural_body(neur_params::NeuronParams;
                              V::Union{Nothing,AbstractMatrix}=nothing,
                              dists::Union{Nothing,AbstractMatrix}=nothing,
                              Rtear::Union{Nothing,AbstractVector}=nothing,
                              rng::AbstractRNG=Random.default_rng())
    pwr    = 1
    nucoff = 3.0
    if V === nothing
        V = spiral_sample_sphere(neur_params.n_samps)
    end
    N = size(V, 1)

    # Teardrop mean (pyramidal cells) or sphere (other types).
    Vtear = if neur_params.neur_type === :pyr
        teardrop_projection(V, 1)
    elseif neur_params.neur_type === :peanut
        teardrop_projection(V, 2)
    else
        copy(V)
    end

    # Covariance based on geodesic (arc-length) distance on the sphere.
    if dists === nothing
        dists = Matrix{Float64}(undef, N, N)
        @inbounds for j in 1:N, i in 1:N
            dx = V[i, 1] - V[j, 1]
            dy = V[i, 2] - V[j, 2]
            dz = V[i, 3] - V[j, 3]
            chord = sqrt(dx^2 + dy^2 + dz^2)
            dists[i, j] = 2 * asin(clamp(chord / 2, -1.0, 1.0))
        end
    end
    C = neur_params.p_scale .* exp.(-(dists ./ neur_params.l_scale) .^ pwr)

    Rt = if Rtear !== nothing
        collect(Rtear)
    elseif neur_params.neur_type === :pyr || neur_params.neur_type === :peanut
        [sqrt(Vtear[i, 1]^2 + Vtear[i, 2]^2 + Vtear[i, 3]^2) for i in 1:N]
    else
        ones(N)
    end

    # Combat ill-conditioning: shift diagonal by |min eigval| if negative.
    min_eig = eigmin(Symmetric((C + C') / 2))
    if min_eig < 0
        C = C + abs(min_eig) * 1.03 * I
    end

    # Sample x_base ~ |N(0, C)|.
    L = cholesky(Symmetric((C + C') / 2)).L
    x_base = abs.(L * randn(rng, N))

    x_bnds = neur_params.exts .* neur_params.avg_rad
    x = x_base .- mean(x_base) .+ neur_params.avg_rad .* Rt
    xmin = min(minimum(x), x_bnds[1])
    xmax = max(maximum(x), x_bnds[2])
    x = (x_bnds[2] - x_bnds[1]) .* (x .- xmin) ./ (xmax - xmin) .+ x_bnds[1]

    x2 = if neur_params.neur_type === :pyr
        xb = x_base .- mean(x_base) .+ neur_params.avg_rad
        xmin = min(minimum(xb), x_bnds[1])
        xmax = max(maximum(xb), x_bnds[2])
        (x_bnds[2] - x_bnds[1]) .* (xb .- xmin) ./ (xmax - xmin) .+ x_bnds[1]
    else
        copy(x)
    end

    # Per-axis eccentricities, normalised to preserve volume.
    eccens = 1 .+ neur_params.eccen .* (rand(rng, 3) .- [0.5, 0.5, 0.0])
    eccens = eccens ./ prod(eccens)^(1 / 3)

    Vetear = if neur_params.neur_type === :pyr
        Vtear .* reshape(eccens, 1, 3)
    else
        V .* reshape(eccens, 1, 3)
    end
    # Normalise so mean-squared radius is unity.
    msr = sqrt(sum(sum(Vetear .^ 2; dims=2)) / size(Vetear, 1))
    Vetear ./= msr

    Vcell = Vetear .* reshape(x, N, 1)
    Vcell[:, 3] .-= nucoff

    Vnorms  = [sqrt(sum(@view(Vcell[i, :]) .^ 2)) for i in 1:N]
    # Nucleus is the soma shape mirrored in z and rescaled.
    Vnuc = V .* reshape([1.0, 1.0, -1.0], 1, 3) .* reshape(x2, N, 1)
    Vnorms2 = [sqrt(sum(@view(Vnuc[i, :]) .^ 2)) for i in 1:N]
    Vnorms2 = neur_params.nexts[2] .* (neur_params.nexts[1] .* (Vnorms2 .- minimum(Vnorms2)) .+
                                         (1 - neur_params.nexts[1]) .* maximum(Vnorms2))
    Vnorms2 = Vnorms2 .+ minimum(Vnorms .- Vnorms2) .- neur_params.min_thic[1]
    Vnuc = Vnuc .* reshape(eccens, 1, 3) .* reshape(Vnorms2 ./
                  [sqrt(sum(@view(Vnuc[i, :]) .^ 2)) for i in 1:N], N, 1)

    lat_ang  = rand(rng) * 2π
    lat_shft = (1 - abs(rand(rng) - rand(rng))) * neur_params.min_thic[2] *
               [sin(lat_ang), cos(lat_ang)]

    Vcell[:, 3] .+= nucoff
    Vnuc[:, 1] .+= lat_shft[1]
    Vnuc[:, 2] .+= lat_shft[2]
    Vnuc[:, 3] .+= nucoff

    # Nucleus volume normalisation (upstream uses convhull volume; we use
    # the radial-mean cube as an analytic proxy that matches to within a
    # few percent for star-shaped meshes).
    if !isempty(neur_params.nuc_rad)
        VnucSz = _star_volume(Vnuc)
        target = (4 / 3) * π * neur_params.nuc_rad[1]^3
        if VnucSz > 0
            factor = (target / VnucSz)^(1 / 3)
            if length(neur_params.nuc_rad) > 1 && neur_params.nuc_rad[2] != 0
                factor = factor^(1 / neur_params.nuc_rad[2])
            end
            Vnuc .*= factor
        end
    end

    # Random rigid rotation (degrees) bounded by `max_ang`.
    max_ang = neur_params.max_ang
    a = -abs(max_ang) .+ 2 * abs(max_ang) .* rand(rng, 3)
    Rx = [1 0 0; 0 cosd(a[1]) -sind(a[1]); 0 sind(a[1]) cosd(a[1])]
    Ry = [cosd(a[2]) 0 sind(a[2]); 0 1 0; -sind(a[2]) 0 cosd(a[2])]
    Rz = [cosd(a[3]) -sind(a[3]) 0; sind(a[3]) cosd(a[3]) 0; 0 0 1]
    R  = Rx * Ry * Rz
    Vcell = Vcell * R
    Vnuc  = Vnuc  * R
    return Vcell, Vnuc, a
end

# Star-shape volume estimate: sum of (1/3) r³ Ω for each vertex, where Ω is
# the per-vertex solid-angle approximation (4π / N for ~uniform sphere
# samplings).
function _star_volume(V::AbstractMatrix{<:Real})
    N = size(V, 1)
    ω = 4π / N
    s = 0.0
    @inbounds for i in 1:N
        r2 = V[i, 1]^2 + V[i, 2]^2 + V[i, 3]^2
        s += r2 * sqrt(r2)   # r^3
    end
    return s * ω / 3
end

# Local `mean` to avoid pulling Statistics into NAOMi's [deps].
function mean(x::AbstractArray)
    s = zero(eltype(x))
    @inbounds for v in x
        s += v
    end
    return s / length(x)
end

# ---------------------------------------------------------------------------
# Point-in-soma test (star-shaped radial check)
# ---------------------------------------------------------------------------

"""
    point_in_soma(Vsurface::AbstractMatrix, p::AbstractVector,
                  center::AbstractVector) -> Bool

Test whether `p` lies inside the star-shaped soma defined by surface
vertices `Vsurface` centred on `center`. Finds the surface vertex whose
direction (relative to `center`) is closest to that of `p`, then
compares the radii.

This replaces upstream's `intriangulation(Vcell, Tri, idx_tri)`. Equivalent
to a nearest-vertex radial test on a star-shaped boundary mesh; accuracy
scales with the angular sample density (`neur_params.n_samps ≈ 200` ⇒
~7° per neighbour).
"""
function point_in_soma(Vsurface::AbstractMatrix{<:Real},
                       p::AbstractVector{<:Real},
                       center::AbstractVector{<:Real})
    dpx = p[1] - center[1]
    dpy = p[2] - center[2]
    dpz = p[3] - center[3]
    rp2 = dpx^2 + dpy^2 + dpz^2
    rp2 == 0 && return true
    rp = sqrt(rp2)
    # find vertex maximizing the cosine of the angle (i.e. closest in
    # angle) between (p - center) and (V - center).
    N = size(Vsurface, 1)
    best_cos = -Inf
    best_r  = Inf
    @inbounds for i in 1:N
        vx = Vsurface[i, 1] - center[1]
        vy = Vsurface[i, 2] - center[2]
        vz = Vsurface[i, 3] - center[3]
        vr = sqrt(vx^2 + vy^2 + vz^2)
        vr == 0 && continue
        c = (dpx * vx + dpy * vy + dpz * vz) / (rp * vr)
        if c > best_cos
            best_cos = c
            best_r  = vr
        end
    end
    return rp <= best_r
end

# ---------------------------------------------------------------------------
# Sample many soma locations + shapes
# ---------------------------------------------------------------------------

"""
    sample_dense_neurons(neur_params, vol_params, neur_ves;
                         rng=Random.default_rng())
        -> (neur_locs, Vcells, Vnucs, rotations)

Sample locations and shapes for up to `vol_params.N_neur` somata,
rejecting positions too close to existing vasculature (`neur_ves`) or
to previously-placed somata (within `vol_params.min_dist`). Returns:

- `neur_locs` — `K × 3` matrix of soma centre positions (microns).
- `Vcells` — length-`K` vector of `N × 3` soma surface meshes (microns,
  shifted to each `neur_locs[k, :]`).
- `Vnucs` — length-`K` vector of `N × 3` nucleus surface meshes.
- `rotations` — `K × 3` matrix of per-cell Euler rotations (degrees).

`K` is at most `N_neur` (may be smaller if the volume is too crowded).
Ports `sampleDenseNeurons.m`.
"""
function sample_dense_neurons(neur_params::NeuronParams,
                              vol_params::VolumeParams,
                              neur_ves::AbstractArray{Bool,3};
                              rng::AbstractRNG=Random.default_rng())
    eta = 1.1
    vres    = vol_params.vres
    vol_sz  = vol_params.vol_sz
    vol_depth = Int(round(vol_params.vol_depth * vres))
    N_neur  = vol_params.N_neur

    # Dilate vasculature by min_dist/2 to keep somata clear of vessels.
    rr = max(1, Int(ceil(vol_params.min_dist / 2)))
    ves_trunc = copy(neur_ves)
    _dilate3d_ball!(ves_trunc, rr)

    # Restrict to the brain-volume slab (skip the cortical-light-path slab).
    H = Int(round(vol_sz[1] * vres))
    W = Int(round(vol_sz[2] * vres))
    D = Int(round(vol_sz[3] * vres))
    idx_good = falses(H, W, D)
    @inbounds for k in 1:D, j in 1:W, i in 1:H
        idx_good[i, j, k] = !ves_trunc[i, j, vol_depth + k]
    end
    idx_bad = copy(idx_good)

    V_master = spiral_sample_sphere(neur_params.n_samps)

    # Pre-compute geodesic distance matrix once and reuse across cells.
    N = size(V_master, 1)
    Dmat = Matrix{Float64}(undef, N, N)
    @inbounds for j in 1:N, i in 1:N
        dx = V_master[i, 1] - V_master[j, 1]
        dy = V_master[i, 2] - V_master[j, 2]
        dz = V_master[i, 3] - V_master[j, 3]
        chord = sqrt(dx^2 + dy^2 + dz^2)
        Dmat[i, j] = 2 * asin(clamp(chord / 2, -1.0, 1.0))
    end

    neur_locs = zeros(Float64, 0, 3)
    Vcells = Vector{Matrix{Float64}}()
    Vnucs  = Vector{Matrix{Float64}}()
    rotations = zeros(Float64, 0, 3)

    while sum(idx_good) > 1 && length(Vcells) < N_neur
        V_tmp, Vnuc_tmp, rot = generate_neural_body(neur_params;
                                                   V=V_master,
                                                   dists=Dmat,
                                                   rng=rng)
        push!(Vcells, V_tmp)
        push!(Vnucs,  Vnuc_tmp)
        rotations = vcat(rotations, reshape(rot, 1, 3))

        # Sample a random allowed voxel position weighted uniformly.
        gpool = idx_good
        tot = sum(gpool)
        if tot == 0
            gpool = idx_bad
            tot = sum(gpool)
            tot == 0 && break
        end
        target = rand(rng, 1:tot)
        chosen = (0, 0, 0)
        running = 0
        @inbounds for k in 1:D, j in 1:W, i in 1:H
            if gpool[i, j, k]
                running += 1
                if running == target
                    chosen = (i, j, k)
                    break
                end
            end
        end
        # Convert voxel index → micron coordinate (linspace 0..vol_sz).
        i, j, k = chosen
        new_pt = if N_neur == 1
            Float64.(vol_sz ./ 2)
        else
            Float64[(i - 1) / (H - 1) * vol_sz[1],
                    (j - 1) / (W - 1) * vol_sz[2],
                    (k - 1) / (D - 1) * vol_sz[3]]
        end
        neur_locs = vcat(neur_locs, reshape(new_pt, 1, 3))

        # Mark a sphere of radius eta*min_dist around new_pt as occupied.
        rsq  = (eta * vol_params.min_dist)^2
        rsq_bad = vol_params.min_dist^2
        @inbounds for kk in 1:D
            zµ = (kk - 1) / (D - 1) * vol_sz[3]
            dz = zµ - new_pt[3]
            (dz^2 > rsq) && continue
            for jj in 1:W
                yµ = (jj - 1) / (W - 1) * vol_sz[2]
                dy = yµ - new_pt[2]
                (dy^2 + dz^2 > rsq) && continue
                for ii in 1:H
                    xµ = (ii - 1) / (H - 1) * vol_sz[1]
                    dx = xµ - new_pt[1]
                    d2 = dx^2 + dy^2 + dz^2
                    if d2 <= rsq
                        idx_good[ii, jj, kk] = false
                    end
                    if d2 <= rsq_bad
                        idx_bad[ii, jj, kk] = false
                    end
                end
            end
        end
        # Apply the upstream `idx_good = ~(idx_good | ~idx_bad)`:
        # equivalent to `idx_good = (~idx_good) & idx_bad`.
        @inbounds for k in 1:D, j in 1:W, i in 1:H
            idx_good[i, j, k] = (!idx_good[i, j, k]) & idx_bad[i, j, k]
        end
    end

    K = length(Vcells)
    # Shift soma + nucleus meshes to their assigned centres.
    Vcells_out = Vector{Matrix{Float64}}(undef, K)
    Vnucs_out  = Vector{Matrix{Float64}}(undef, K)
    for k in 1:K
        loc = neur_locs[k, :]
        Vcells_out[k] = Vcells[k] .+ reshape(loc, 1, 3)
        Vnucs_out[k]  = Vnucs[k]  .+ reshape(loc, 1, 3)
    end

    return neur_locs, Vcells_out, Vnucs_out, rotations
end

# A 3-D dilation by a ball; used to fence somas off from vasculature.
function _dilate3d_ball!(vol::AbstractArray{Bool,3}, r::Integer)
    src = copy(vol)
    fill!(vol, false)
    H, W, D = size(vol)
    r2 = Float64(r)^2
    @inbounds for kk in 1:D, jj in 1:W, ii in 1:H
        src[ii, jj, kk] || continue
        for dk in -r:r, dj in -r:r, di in -r:r
            (di^2 + dj^2 + dk^2) <= r2 || continue
            i = ii + di
            j = jj + dj
            k = kk + dk
            (1 <= i <= H && 1 <= j <= W && 1 <= k <= D) || continue
            vol[i, j, k] = true
        end
    end
    return vol
end

# ---------------------------------------------------------------------------
# Rasterize soma/nucleus meshes into voxel volumes
# ---------------------------------------------------------------------------

"""
    generate_neural_volume(neur_params, vol_params, neur_locs, Vcells, Vnucs,
                           neur_ves) -> (neur_soma, neur_vol, gp_nuc, gp_soma)

Rasterize each cell's surface mesh into a `vol_sz × vres` voxel volume.
Returns:

- `neur_soma::Array{UInt16,3}` — voxel-wise neuron index (0 = empty).
- `neur_vol::Array{Float32,3}` — base fluorescence map (initialised with
  `neur_params.nuc_fluorsc` inside nuclei).
- `gp_nuc::Vector{Tuple{Vector{Int32},Float64}}` — per-cell `(indices, val)`
  pair for nucleus voxels.
- `gp_soma::Vector{Vector{Int32}}` — per-cell linear-index list of soma
  voxels (excluding nucleus).

Replaces upstream `intriangulation(Vcell, Tri, ...)` with a star-shape
radial test (`point_in_soma`). Ports `generateNeuralVolume.m`.
"""
function generate_neural_volume(neur_params::NeuronParams,
                                vol_params::VolumeParams,
                                neur_locs::AbstractMatrix{<:Real},
                                Vcells::AbstractVector,
                                Vnucs::AbstractVector,
                                neur_ves::AbstractArray{Bool,3})
    vres = vol_params.vres
    vol_sz = vol_params.vol_sz
    H = Int(round(vol_sz[1] * vres))
    W = Int(round(vol_sz[2] * vres))
    D = Int(round(vol_sz[3] * vres))
    vol_depth = Int(round(vol_params.vol_depth * vres))

    neur_soma = zeros(UInt16, H, W, D)
    neur_vol  = zeros(Float32, H, W, D)
    K = length(Vcells)
    gp_nuc  = Vector{Tuple{Vector{Int32}, Float64}}(undef, K)
    gp_soma = Vector{Vector{Int32}}(undef, K)
    taken = falses(H, W, D)
    @inbounds for k in 1:D, j in 1:W, i in 1:H
        taken[i, j, k] = neur_ves[i, j, vol_depth + k]
    end

    linsub(i, j, k) = Int32(LinearIndices((H, W, D))[i, j, k])

    for kk in 1:K
        center = neur_locs[kk, :]
        # Compute soma bbox (in voxel space).
        Vc = Vcells[kk]
        Vn = Vnucs[kk]
        max_ext = 0.0
        for i in axes(Vc, 1)
            r = sqrt(sum((@view(Vc[i, :]) .- center) .^ 2))
            r > max_ext && (max_ext = r)
        end
        mER = Int(ceil(max_ext * vres))
        idxC = (Int(round(vres * center[1])),
                Int(round(vres * center[2])),
                Int(round(vres * center[3])))
        ix_lo = max(1, idxC[1] - mER); ix_hi = min(idxC[1] + mER, H)
        iy_lo = max(1, idxC[2] - mER); iy_hi = min(idxC[2] + mER, W)
        iz_lo = max(1, idxC[3] - mER); iz_hi = min(idxC[3] + mER, D)

        soma_indices = Int32[]
        nuc_indices  = Int32[]
        sizehint!(soma_indices, (ix_hi - ix_lo + 1) * (iy_hi - iy_lo + 1))

        @inbounds for kz in iz_lo:iz_hi, jy in iy_lo:iy_hi, ix in ix_lo:ix_hi
            taken[ix, jy, kz] && continue
            # voxel center coords (in microns):
            p = (ix / vres + 0.5 / vres,
                 jy / vres + 0.5 / vres,
                 kz / vres + 0.5 / vres)
            d2 = (p[1] - center[1])^2 + (p[2] - center[2])^2 + (p[3] - center[3])^2
            d2 > max_ext^2 && continue
            in_cell = point_in_soma(Vc, [p[1], p[2], p[3]], center)
            in_cell || continue
            in_nuc  = point_in_soma(Vn, [p[1], p[2], p[3]], center)
            li = linsub(ix, jy, kz)
            if in_nuc
                push!(nuc_indices, li)
                neur_vol[ix, jy, kz] = Float32(neur_params.nuc_fluorsc)
            else
                push!(soma_indices, li)
                neur_soma[ix, jy, kz] = UInt16(kk)
            end
            taken[ix, jy, kz] = true
        end
        gp_nuc[kk]  = (nuc_indices, Float64(neur_params.nuc_fluorsc))
        gp_soma[kk] = soma_indices
    end
    return neur_soma, neur_vol, gp_nuc, gp_soma
end

# ---------------------------------------------------------------------------
# Filter visible somata (z within PSF + soma radius of the imaging plane)
# ---------------------------------------------------------------------------

"""
    isolate_visible_somas(locs, psf, vol_params, neur_params; thresh=0)
        -> Vector{Int}

Return the indices of somata whose centre `z` is within
`thresh + psfHalfWidth/2 + neur_params.avg_rad` of the imaging plane
(volume midplane). `psf` must be a 3-D array; the half-width is taken
from [`width_estimate_3d`](@ref) (Chunk 7). Ports `isolateVisibleSomas.m`.
"""
function isolate_visible_somas(locs::AbstractMatrix{<:Real},
                               psf::AbstractArray{<:Real,3},
                               vol_params::VolumeParams,
                               neur_params::NeuronParams;
                               thresh::Real=0)
    vMid = vol_params.vol_sz[3] / 2
    zLocs = locs[:, end]
    psf_half = width_estimate_3d(psf)[end]
    t = thresh + psf_half / 2 + neur_params.avg_rad
    return findall(abs.(zLocs .- vMid) .< t)
end

# ---------------------------------------------------------------------------
# 3-D Gaussian process via FFT (small helper used by the fluorescence chunk)
# ---------------------------------------------------------------------------

"""
    masked_3d_gp(grid_sz, l_scale, p_scale, mu;
                 bin_mask=1, l_weights=1, rng=Random.default_rng())
        -> Array{Float32,3}

Draw a sample from a 3-D Gaussian process with mean `mu`, length scale
`l_scale` (scalar or 3-vector), variance `p_scale`, on a grid of size
`grid_sz`. The kernel is applied in the frequency domain (FFT-based
sampling), which is fast for large grids. Optional `bin_mask` is
applied multiplicatively at the end. Ports `masked_3DGP_v2.m`.
"""
function masked_3d_gp(grid_sz, l_scale, p_scale::Real, mu;
                      bin_mask=1.0f0, l_weights=1.0,
                      rng::AbstractRNG=Random.default_rng())
    sz = length(grid_sz) == 1 ?
        (Int(grid_sz), Int(grid_sz), Int(grid_sz)) :
        (Int(grid_sz[1]), Int(grid_sz[2]), Int(grid_sz[3]))
    L = if l_scale isa AbstractMatrix
        Float64.(l_scale)
    elseif l_scale isa AbstractVector
        reshape(Float64.(l_scale), 1, length(l_scale))
    else
        reshape(Float64[l_scale l_scale l_scale], 1, 3)
    end
    if size(L, 2) == 1
        L = hcat(L, L, L)
    end
    nl = size(L, 1)
    lw = l_weights isa AbstractArray ? Float64.(l_weights) :
         fill(Float64(l_weights), nl)

    wmx = π / 2
    gx = reshape(range(-wmx, wmx; length=sz[1]) .^ 2, :, 1, 1)
    gy = reshape(range(-wmx, wmx; length=sz[2]) .^ 2, 1, :, 1)
    gz = reshape(range(-wmx, wmx; length=sz[3]) .^ 2, 1, 1, :)
    gp = zeros(ComplexF32, sz)
    for i in 1:nl
        ker = (exp.(-gx .* L[i, 1]^2) .* exp.(-gy .* L[i, 2]^2)) .* exp.(-gz .* L[i, 3]^2)
        TMP = ComplexF32.(randn(rng, Float32, sz) .+ 1im .* randn(rng, Float32, sz))
        scale = Float32(lw[i] * sqrt(L[i, 1] * L[i, 2] * L[i, 3]))
        gp .+= scale .* TMP .* ker
    end
    gp_real = sqrt(Float32(prod(sz))) .* real.(ifftshift(ifft(ifftshift(gp))))
    norm = Float32(p_scale * (2^4.5 / π^1.5) / sqrt(nl))
    return norm .* bin_mask .* gp_real .+ Float32(mu)
end
