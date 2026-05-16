# Dendrite generation.
#
# Ported from the upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - dendrite_dijkstra2.m + dendrite_dijkstra_cpp.cpp → dendrite_dijkstra_grid
#   - getDendritePath2.m                               → get_dendrite_path
#   - dendrite_randomwalk2.m + ..._cpp.cpp             → dendrite_random_walk
#   - dilateDendritePathAll.m                          → dilate_dendrite_paths_all
#   - growNeuronDendrites.m                            → grow_neuron_dendrites!
#   - growApicalDendrites.m                            → grow_apical_dendrites!
#   - smoothCellBody.m                                 → smooth_cell_body
#   - setCellFluoresence.m                             → set_cell_fluorescence

export dendrite_dijkstra_grid, get_dendrite_path, dendrite_random_walk,
       dilate_dendrite_paths_all, grow_neuron_dendrites!,
       grow_apical_dendrites!, smooth_cell_body, set_cell_fluorescence

# ---------------------------------------------------------------------------
# Binary min-heap (small, hand-rolled — avoids adding `DataStructures`)
# ---------------------------------------------------------------------------

mutable struct _MinHeap
    keys::Vector{Float32}
    vals::Vector{Int}
end

_MinHeap() = _MinHeap(Float32[], Int[])
Base.isempty(h::_MinHeap) = isempty(h.keys)

function _push!(h::_MinHeap, k::Float32, v::Int)
    push!(h.keys, k)
    push!(h.vals, v)
    i = length(h.keys)
    @inbounds while i > 1
        p = i >> 1
        if h.keys[p] > h.keys[i]
            h.keys[p], h.keys[i] = h.keys[i], h.keys[p]
            h.vals[p], h.vals[i] = h.vals[i], h.vals[p]
            i = p
        else
            break
        end
    end
end

function _pop!(h::_MinHeap)
    n = length(h.keys)
    n == 0 && throw(BoundsError(h))
    k0 = h.keys[1]
    v0 = h.vals[1]
    if n == 1
        pop!(h.keys); pop!(h.vals)
        return k0, v0
    end
    h.keys[1] = pop!(h.keys)
    h.vals[1] = pop!(h.vals)
    nn = length(h.keys)
    i = 1
    @inbounds while true
        l = 2 * i
        r = l + 1
        smallest = i
        if l <= nn && h.keys[l] < h.keys[smallest]
            smallest = l
        end
        if r <= nn && h.keys[r] < h.keys[smallest]
            smallest = r
        end
        smallest == i && break
        h.keys[i], h.keys[smallest] = h.keys[smallest], h.keys[i]
        h.vals[i], h.vals[smallest] = h.vals[smallest], h.vals[i]
        i = smallest
    end
    return k0, v0
end

# ---------------------------------------------------------------------------
# Grid Dijkstra on a 6-directional edge-weight array
# ---------------------------------------------------------------------------

# Core grid-Dijkstra with reusable buffers and optional early termination.
#
# `distance` / `pathfrom` are reset only over the `touched` voxels left by
# the previous call, then reused — the caller allocates them once and
# threads them through every per-cell call. `heap` is likewise reused (its
# capacity is retained across calls).
#
# With a non-empty `targets` list (linear voxel indices) the search stops
# as soon as every target voxel has been settled. This is exact: Dijkstra
# settles voxels in nondecreasing distance order, so once a voxel is
# popped its `distance` and `pathfrom` entry are already final — and
# `get_dendrite_path` only ever walks back from target voxels.
function _dijkstra!(distance::AbstractArray{Float32,3},
                    pathfrom::AbstractArray{Int,3},
                    heap::_MinHeap, touched::Vector{Int},
                    M::AbstractArray{<:Real,4},
                    root::NTuple{3,<:Integer},
                    targets::AbstractVector{<:Integer})
    H, W, D = size(M, 1), size(M, 2), size(M, 3)
    size(M, 4) == 6 || throw(ArgumentError("M must have size (H, W, D, 6)"))
    lin = LinearIndices((H, W, D))
    cart = CartesianIndices((H, W, D))

    @inbounds for li in touched
        distance[li] = Inf32
        pathfrom[li] = 0
    end
    empty!(touched)
    empty!(heap.keys)
    empty!(heap.vals)

    root_li = lin[root...]
    distance[root_li] = 0.0f0
    pathfrom[root_li] = 0
    push!(touched, root_li)
    _push!(heap, 0.0f0, root_li)

    early = !isempty(targets)
    target_set = Set{Int}()
    for t in targets
        push!(target_set, Int(t))
    end
    n_left = length(target_set)

    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
               (0, -1, 0), (0, 0, 1), (0, 0, -1))

    while !isempty(heap)
        cn_d, cn_idx = _pop!(heap)
        cn_d > distance[cn_idx] && continue
        if early && (cn_idx in target_set)
            delete!(target_set, cn_idx)
            n_left -= 1
            n_left == 0 && break
        end
        c = cart[cn_idx]
        @inbounds for d in 1:6
            off = offsets[d]
            ni = c[1] + off[1]
            nj = c[2] + off[2]
            nk = c[3] + off[3]
            (1 <= ni <= H && 1 <= nj <= W && 1 <= nk <= D) || continue
            w = M[ni, nj, nk, d]
            isfinite(w) || continue
            ndist = cn_d + Float32(w)
            nidx = lin[ni, nj, nk]
            if ndist < distance[nidx]
                distance[nidx] == Inf32 && push!(touched, nidx)
                distance[nidx] = ndist
                pathfrom[nidx] = cn_idx
                _push!(heap, ndist, nidx)
            end
        end
    end
    return distance, pathfrom
end

"""
    dendrite_dijkstra_grid(M::AbstractArray{<:Real,4}, root::NTuple{3,<:Integer})
        -> (distance, pathfrom)

Single-source shortest path on a 3-D voxel grid where `M[i, j, k, d]` is
the cost of *entering* voxel `(i, j, k)` from direction `d ∈ 1:6`
(directions: `+x, -x, +y, -y, +z, -z`).

Returns:

- `distance::Array{Float32,3}` — shortest path length from `root` to each
  voxel (`Inf` if unreachable).
- `pathfrom::Array{Int,3}` — linear index (1-based) of the predecessor
  voxel on the shortest path (`0` for `root` and unreachable voxels).

Ports `dendrite_dijkstra2.m` + the `dendrite_dijkstra_cpp.cpp` MEX
kernel; uses a hand-rolled binary min-heap. The dendrite-growth callers
use an internal early-terminating variant that stops once the sampled
endpoint voxels are settled; this exported entry point runs the full
single-source search.
"""
function dendrite_dijkstra_grid(M::AbstractArray{<:Real,4},
                                root::NTuple{3,<:Integer})
    H, W, D = size(M, 1), size(M, 2), size(M, 3)
    size(M, 4) == 6 || throw(ArgumentError("M must have size (H, W, D, 6)"))
    distance = fill(Float32(Inf), H, W, D)
    pathfrom = zeros(Int, H, W, D)
    _dijkstra!(distance, pathfrom, _MinHeap(), Int[], M, root, Int[])
    return distance, pathfrom
end

# ---------------------------------------------------------------------------
# Walk the predecessor map from an endpoint back to root
# ---------------------------------------------------------------------------

"""
    get_dendrite_path(pathfrom::AbstractArray{<:Integer,3},
                      node::NTuple{3,<:Integer},
                      root::NTuple{3,<:Integer})
        -> Vector{NTuple{3,Int}}

Walk back from `node` to `root` through the 3-D `pathfrom` predecessor
map (linear indices). Returns an empty vector if `node` is unreachable.
Ports `getDendritePath2.m`.
"""
function get_dendrite_path(pathfrom::AbstractArray{<:Integer,3},
                           node::NTuple{3,<:Integer},
                           root::NTuple{3,<:Integer})
    H, W, D = size(pathfrom)
    cart = CartesianIndices((H, W, D))
    path = NTuple{3,Int}[]
    cur = (Int(node[1]), Int(node[2]), Int(node[3]))
    push!(path, cur)
    safety = H * W * D + 1
    while cur != root && safety > 0
        safety -= 1
        idx = pathfrom[cur...]
        idx == 0 && return NTuple{3,Int}[]
        c = cart[idx]
        cur = (c[1], c[2], c[3])
        push!(path, cur)
    end
    return path
end

# ---------------------------------------------------------------------------
# Greedy directed random walk
# ---------------------------------------------------------------------------

"""
    dendrite_random_walk(M::AbstractArray{<:Real,3}, root, ends;
                         distsc=1.0, maxlength=2000, fillweight=1.0,
                         maxel=8, minlength=10)
        -> Vector{NTuple{3,Int}}

Greedy locally-best directed walk from `root` toward `ends` in an
occupancy field `M`. Mutates `M` along the accepted path
(`M[p] += fillweight`). Returns the path; empty if `minlength` was not
reached. Ports `dendrite_randomwalk2.m` + `dendrite_randomwalk_cpp.cpp`.
"""
function dendrite_random_walk(M::AbstractArray{<:Real,3},
                              root::NTuple{3,<:Integer},
                              ends::NTuple{3,<:Integer};
                              distsc::Real=1.0,
                              maxlength::Integer=2000,
                              fillweight::Real=1.0,
                              maxel::Integer=8,
                              minlength::Integer=10)
    H, W, D = size(M)
    cur = (Int(root[1]), Int(root[2]), Int(root[3]))
    target = (Int(ends[1]), Int(ends[2]), Int(ends[3]))
    maxfill = Float32(maxel * fillweight)
    path = NTuple{3,Int}[]
    Mvals = Float32[]
    bglength = 0
    for i in 1:maxlength
        dx = cur[1] - target[1]
        dy = cur[2] - target[2]
        dz = cur[3] - target[3]
        dist = max(1e-6, sqrt(dx^2 + dy^2 + dz^2) / distsc)
        dxn = dx / dist
        dyn = dy / dist
        dzn = dz / dist

        jmin = 0
        minmat = Float32(Inf)
        if cur[1] != H
            v = Float32(M[cur[1] + 1, cur[2], cur[3]] + dxn)
            v < minmat && (minmat = v; jmin = 1)
        end
        if cur[2] != W
            v = Float32(M[cur[1], cur[2] + 1, cur[3]] + dyn)
            v < minmat && (minmat = v; jmin = 2)
        end
        if cur[3] != D
            v = Float32(M[cur[1], cur[2], cur[3] + 1] + dzn)
            v < minmat && (minmat = v; jmin = 3)
        end
        if cur[1] != 1
            v = Float32(M[cur[1] - 1, cur[2], cur[3]] - dxn)
            v < minmat && (minmat = v; jmin = 4)
        end
        if cur[2] != 1
            v = Float32(M[cur[1], cur[2] - 1, cur[3]] - dyn)
            v < minmat && (minmat = v; jmin = 5)
        end
        if cur[3] != 1
            v = Float32(M[cur[1], cur[2], cur[3] - 1] - dzn)
            v < minmat && (minmat = v; jmin = 6)
        end

        if minmat < maxfill
            nx = jmin == 1 ? cur[1] + 1 : jmin == 4 ? cur[1] - 1 : cur[1]
            ny = jmin == 2 ? cur[2] + 1 : jmin == 5 ? cur[2] - 1 : cur[2]
            nz = jmin == 3 ? cur[3] + 1 : jmin == 6 ? cur[3] - 1 : cur[3]
            cur = (nx, ny, nz)
            push!(path, cur)
            push!(Mvals, minmat)
            M[cur...] = Float32(Inf)
            (cur[1] == 1 || cur[2] == 1 || cur[3] == 1 ||
             cur[1] == H || cur[2] == W || cur[3] == D) && break
            cur == target && break
        else
            bglength = i - 1
            break
        end
    end
    bglength = bglength <= 0 ? length(path) : bglength
    if bglength >= minlength
        for i in 1:bglength
            if Mvals[i] < maxfill
                M[path[i]...] = Mvals[i] + Float32(fillweight)
            end
        end
        return path[1:bglength]
    else
        for i in 1:bglength
            M[path[i]...] = Mvals[i]
        end
        return NTuple{3,Int}[]
    end
end

# ---------------------------------------------------------------------------
# Iterative path dilation
# ---------------------------------------------------------------------------

"""
    dilate_dendrite_paths_all(paths, pathnums, obstruction;
                              maxDist::Integer=20, rng=Random.default_rng())
        -> (paths_out, pathnums_out)

Grow each labelled dendrite path (encoded with values `> 1` in `paths`)
outward into adjacent empty voxels (`paths == 0`), preserving labels in
`pathnums`. Respects an `obstruction` map. Ports `dilateDendritePathAll.m`.
"""
function dilate_dendrite_paths_all(paths::AbstractArray{<:Real,3},
                                   pathnums::AbstractArray{<:Integer,3},
                                   obstruction::AbstractArray{<:Real,3};
                                   maxDist::Integer=20,
                                   rng::AbstractRNG=Random.default_rng())
    H, W, D = size(paths)
    paths_out = Float32.(paths)
    pathnums_out = Int.(pathnums)
    for i in eachindex(paths_out)
        if obstruction[i] != 0
            paths_out[i] = NaN
        end
    end

    lin = LinearIndices((H, W, D))
    cart = CartesianIndices((H, W, D))

    offsets = NTuple{3,Int}[]
    dist2 = Int[]
    for dz in -maxDist:maxDist, dy in -maxDist:maxDist, dx in -maxDist:maxDist
        d2 = dx^2 + dy^2 + dz^2
        d2 == 0 && continue
        d2 <= maxDist^2 || continue
        push!(offsets, (dx, dy, dz))
        push!(dist2, d2)
    end
    perm = sortperm(dist2)
    offsets = offsets[perm]
    dist2 = dist2[perm]
    shells = Vector{UnitRange{Int}}()
    i = 1
    while i <= length(dist2)
        j = i
        while j <= length(dist2) && dist2[j] == dist2[i]
            j += 1
        end
        push!(shells, i:(j - 1))
        i = j
    end
    neigh6 = ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
              (0, -1, 0), (0, 0, 1), (0, 0, -1))

    for s in 1:length(shells)
        any_remaining = false
        for c in eachindex(paths_out)
            if isfinite(paths_out[c]) && paths_out[c] > 1
                any_remaining = true
                break
            end
        end
        any_remaining || break
        hot = Int[]
        for ix in eachindex(paths_out)
            if isfinite(paths_out[ix]) && paths_out[ix] > 1
                push!(hot, ix)
            end
        end
        shell_offs = offsets[shells[s]]
        for j in hot
            paths_out[j] > 1 || continue
            jc = cart[j]
            candidates = Int[]
            for off in shell_offs
                ii = jc[1] + off[1]
                jj = jc[2] + off[2]
                kk = jc[3] + off[3]
                (1 <= ii <= H && 1 <= jj <= W && 1 <= kk <= D) || continue
                isfinite(paths_out[ii, jj, kk]) || continue
                paths_out[ii, jj, kk] == 0 || continue
                same_label = false
                for nb in neigh6
                    ai = ii + nb[1]; aj = jj + nb[2]; ak = kk + nb[3]
                    (1 <= ai <= H && 1 <= aj <= W && 1 <= ak <= D) || continue
                    if pathnums_out[ai, aj, ak] == pathnums_out[j]
                        same_label = true
                        break
                    end
                end
                same_label && push!(candidates, lin[ii, jj, kk])
            end
            while paths_out[j] > 1 && !isempty(candidates)
                r = max(1, ceil(Int, rand(rng) * length(candidates)))
                pidx = candidates[r]
                deleteat!(candidates, r)
                paths_out[j] -= 1
                paths_out[pidx] = 1
                pathnums_out[pidx] = pathnums_out[j]
            end
        end
    end

    for i in eachindex(paths_out)
        if !isfinite(paths_out[i])
            paths_out[i] = 0
        end
    end
    return paths_out, pathnums_out
end

# ---------------------------------------------------------------------------
# Build the 6-direction edge-weight array for Dijkstra-on-grid
# ---------------------------------------------------------------------------

# M[i,j,k,d] = cost to enter voxel (i,j,k) from direction d.
function _build_dendrite_M(dims::NTuple{3,Int},
                          obstruction::AbstractArray{<:Real,3},
                          dweight::Real, bweight::Real;
                          rng::AbstractRNG=Random.default_rng())
    H, W, D = dims
    M = ones(Float32, H, W, D, 6) .+ Float32(dweight) .* rand(rng, Float32, H, W, D, 6)
    M[1,   :, :, 1] .= Inf32
    M[end, :, :, 2] .= Inf32
    M[:, 1,   :, 3] .= Inf32
    M[:, end, :, 4] .= Inf32
    M[:, :, 1,   5] .= Inf32
    M[:, :, end, 6] .= Inf32
    fillfrac = Float32.(obstruction .> 0)
    pen = -Float32(bweight) .* log.(max.(1f-6, 1 .- 2 .* max.(0.0f0, fillfrac .- 0.5f0)))
    for d in 1:6
        M[:, :, :, d] .+= pen
    end
    return M
end

# ---------------------------------------------------------------------------
# Per-neuron dendrite growth (basal tree + local apical tree)
# ---------------------------------------------------------------------------

"""
    grow_neuron_dendrites!(vol_params, dend_params, neur_soma, neur_ves,
                           neur_locs, gp_nuc, gp_soma; rng=Random.default_rng())
        -> (neur_num, dendnum_AD, gp_soma_aug)

Grow basal-tree and local-apical-tree dendrites off each soma. Returns:

- `neur_num::Array{UInt16,3}` — voxel-wise neuron index including
  dendrites.
- `dendnum_AD::Array{UInt16,3}` — voxel-wise apical-dendrite index.
- `gp_soma_aug` — per-cell `(soma=..., smoothed=...)` index pairs.

**Single-stage Dijkstra at fine resolution** (no coarse-then-fine
refinement). Ports `growNeuronDendrites.m`.
"""
function grow_neuron_dendrites!(vol_params::VolumeParams,
                                dend_params::DendriteParams,
                                neur_soma::AbstractArray{<:Integer,3},
                                neur_ves::AbstractArray{Bool,3},
                                neur_locs::AbstractMatrix{<:Real},
                                gp_nuc::AbstractVector,
                                gp_soma::AbstractVector;
                                rng::AbstractRNG=Random.default_rng())
    dtParams       = collect(Float64, dend_params.dtParams)
    atParams       = collect(Float64, dend_params.atParams)
    dweight        = dend_params.dweight
    bweight        = dend_params.bweight
    thicknessScale = dend_params.thicknessScale
    rallexp        = dend_params.rallexp
    vres           = vol_params.vres
    N_neur         = vol_params.N_neur
    vol_sz         = vol_params.vol_sz
    fulldims = (Int(round(vol_sz[1] * vres)),
                Int(round(vol_sz[2] * vres)),
                Int(round(vol_sz[3] * vres)))
    dtParams[2:3] .*= vres
    atParams[2:4] .*= vres
    thicknessScale *= vres * vres

    vol_depth = Int(round(vol_params.vol_depth * vres))
    cellVolume = Float32.(neur_soma)
    obstruction_marker = Float32(vol_params.N_neur + Int(round(vol_params.N_den)) +
                                 vol_params.N_bg + 1)
    for k in 1:fulldims[3], j in 1:fulldims[2], i in 1:fulldims[1]
        if neur_ves[i, j, vol_depth + k]
            cellVolume[i, j, k] = obstruction_marker
        end
    end
    for kk in 1:length(gp_nuc)
        for x in gp_nuc[kk][1]
            cellVolume[x] = Float32(kk)
        end
    end

    neur_num = UInt16.(clamp.(cellVolume, 0, typemax(UInt16)))
    cellVolumeIdx = zeros(Float32, fulldims)
    cellVolumeVal = zeros(Float32, fulldims)
    cellVolumeAD  = falses(fulldims)
    allroots = Int.(ceil.(max.(vres .* neur_locs, 1e-4)))

    dendVar = 0.25
    gp_soma_aug = Vector{Dict{Symbol,Any}}(undef, N_neur)
    for kk in 1:N_neur
        gp_soma_aug[kk] = Dict{Symbol,Any}(:soma => gp_soma[kk], :smoothed => Int32[])
    end
    lin_full = LinearIndices(fulldims)

    # Reusable Dijkstra buffers (see `_dijkstra!`); allocated once.
    dij_dist = fill(Inf32, fulldims)
    dij_path = zeros(Int, fulldims)
    dij_heap = _MinHeap()
    dij_touched = Int[]

    for j in 1:N_neur
        rootL = (clamp(allroots[j, 1], 1, fulldims[1]),
                 clamp(allroots[j, 2], 1, fulldims[2]),
                 clamp(allroots[j, 3], 1, fulldims[3]))
        obstruction = copy(cellVolume)
        for x in gp_soma[j]
            obstruction[x] = 0
        end
        cellBody = Int[Int(x) for x in gp_soma[j]]

        M = _build_dendrite_M(fulldims, obstruction, dweight, bweight; rng=rng)
        # Anchor: connect rootL → soma "corner" (the lowest-linear-index
        # soma voxel) with a 3-D staircase of zero-cost edges, matching
        # upstream `aproot` semantics in `growNeuronDendrites.m`
        # (lines 191-204).
        if !isempty(gp_soma[j])
            ap_li = Int(minimum(gp_soma[j]))
            ap = Tuple(CartesianIndices(fulldims)[ap_li])
            # Walk x from rootL[1] → ap[1] with y=rootL[2], z=rootL[3]
            if ap[1] > rootL[1]
                for v in rootL[1]:ap[1]
                    M[v, rootL[2], rootL[3], 1] = 0f0
                end
            elseif ap[1] < rootL[1]
                for v in ap[1]:rootL[1]
                    M[v, rootL[2], rootL[3], 2] = 0f0
                end
            end
            # Walk y with x=ap[1], z=rootL[3]
            if ap[2] > rootL[2]
                for v in rootL[2]:ap[2]
                    M[ap[1], v, rootL[3], 3] = 0f0
                end
            elseif ap[2] < rootL[2]
                for v in ap[2]:rootL[2]
                    M[ap[1], v, rootL[3], 4] = 0f0
                end
            end
            # Walk z with x=ap[1], y=ap[2]
            if ap[3] > rootL[3]
                for v in rootL[3]:ap[3]
                    M[ap[1], ap[2], v, 5] = 0f0
                end
            elseif ap[3] < rootL[3]
                for v in ap[3]:rootL[3]
                    M[ap[1], ap[2], v, 6] = 0f0
                end
            end
        end
        # Basal-tree endpoints.
        numdt = max(1, Int(round(dtParams[1] + dtParams[5] * randn(rng))))
        endsT = NTuple{3,Int}[]
        for _ in 1:numdt
            flag = true
            distSC = 1.0
            numit = 0
            while flag && numit < 100
                θ = rand(rng) * 2π
                r = sqrt(rand(rng)) * dtParams[2] * distSC
                z = clamp(Int(floor(2 * dtParams[3] * (rand(rng) - 0.5) + rootL[3])),
                          1, fulldims[3])
                x = clamp(Int(floor(r * cos(θ) + rootL[1])), 1, fulldims[1])
                y = clamp(Int(floor(r * sin(θ) + rootL[2])), 1, fulldims[2])
                if obstruction[x, y, z] == 0
                    push!(endsT, (x, y, z))
                    flag = false
                end
                distSC *= 1.01
                numit += 1
            end
            flag && push!(endsT, rootL)
        end
        # Apical-tree endpoints (local).
        nat = Int(atParams[1])
        rootA = (clamp(Int(floor(rootL[1] + 2 * atParams[4] * (rand(rng) - 0.5))),
                       1, fulldims[1]),
                 clamp(Int(floor(rootL[2] + 2 * atParams[4] * (rand(rng) - 0.5))),
                       1, fulldims[2]),
                 clamp(Int(1 + atParams[3]), 1, fulldims[3]))
        endsA = NTuple{3,Int}[]
        for _ in 1:nat
            flag = true
            distSC = 1.0
            numit = 0
            while flag && numit < 100
                θ = rand(rng) * 2π
                r = sqrt(rand(rng)) * atParams[2] * distSC
                z = clamp(Int(floor(2 * atParams[3] * (rand(rng) - 0.5) + rootA[3])),
                          1, fulldims[3])
                x = clamp(Int(floor(r * cos(θ) + rootA[1])), 1, fulldims[1])
                y = clamp(Int(floor(r * sin(θ) + rootA[2])), 1, fulldims[2])
                if obstruction[x, y, z] == 0
                    push!(endsA, (x, y, z))
                    flag = false
                end
                distSC *= 1.01
                numit += 1
            end
            flag && push!(endsA, rootA)
        end
        all_ends = vcat(endsT, endsA)

        # Endpoint sampling above draws no Dijkstra output, so running it
        # before the search leaves the RNG stream untouched. Settling just
        # those endpoints lets the search stop early.
        targets = Int[lin_full[e...] for e in all_ends]
        _dijkstra!(dij_dist, dij_path, dij_heap, dij_touched, M, rootL, targets)
        pathfrom = dij_path

        all_paths = Vector{Vector{NTuple{3,Int}}}(undef, length(all_ends))

        finepathsIdx = zeros(Float32, fulldims)
        finepathsVal = zeros(Float32, fulldims)
        finepathsAD  = falses(fulldims)
        fineIdxs_basal = Set{Int}()
        fineIdxs_apic  = Set{Int}()
        for (ei, e) in enumerate(all_ends)
            path = get_dendrite_path(pathfrom, e, rootL)
            all_paths[ei] = path
            isempty(path) && continue
            dendSz = max(0.0, 1 + dendVar * randn(rng))
            n = length(path)
            pw = ones(Float32, n)
            if n > 2
                for i in 2:(n - 1)
                    p, c, q = path[i - 1], path[i], path[i + 1]
                    d2 = abs(2c[1] - p[1] - q[1]) +
                         abs(2c[2] - p[2] - q[2]) +
                         abs(2c[3] - p[3] - q[3])
                    pw[i] = Float32(dendSz * (1 - (1 - 1 / sqrt(2)) * d2 / 2))
                end
                pw[1] = pw[2]; pw[end] = pw[end - 1]
                pw .= max.(pw, 0f0)
            else
                pw .= Float32(dendSz)
            end
            apical = ei > numdt
            for (i, p) in enumerate(path)
                li = lin_full[p...]
                if apical
                    finepathsVal[li] += pw[i]
                    finepathsAD[li]   = true
                    push!(fineIdxs_apic, li)
                else
                    finepathsVal[li] += pw[i]
                    push!(fineIdxs_basal, li)
                end
            end
        end

        for li in fineIdxs_basal
            v = finepathsVal[li]
            v > 0 && (finepathsVal[li] = Float32(thicknessScale * dtParams[4] *
                                                 v^(1 / rallexp)))
        end
        for li in fineIdxs_apic
            v = finepathsVal[li]
            v > 0 && (finepathsVal[li] = Float32(thicknessScale * atParams[5] *
                                                 v^(1 / rallexp)))
        end

        fineIdxs3 = collect(union(fineIdxs_basal, fineIdxs_apic))
        for li in fineIdxs3
            finepathsIdx[li] = Float32(j)
        end

        cellBodyVec = Int32.(cellBody)
        if !all(isempty, all_paths)
            smoothed = smooth_cell_body(all_paths, cellBodyVec, fulldims)
            cb_set = Set(cellBody)
            for li in smoothed
                if !(Int(li) in cb_set)
                    finepathsIdx[Int(li)] = Float32(j)
                    finepathsVal[Int(li)] += 1f0
                    push!(fineIdxs3, Int(li))
                end
            end
            gp_soma_aug[j][:smoothed] = collect(smoothed)
        end

        for li in cellBody
            finepathsIdx[li] = 0
            finepathsVal[li] = 0
            finepathsAD[li]  = false
        end

        for li in fineIdxs3
            cellVolume[li]    += finepathsIdx[li]
            cellVolumeIdx[li] += finepathsIdx[li]
            cellVolumeVal[li] += finepathsVal[li]
            cellVolumeAD[li]  |= finepathsAD[li]
        end
    end

    # Stochastic rounding → UInt16 thickness.
    thick = zeros(UInt16, fulldims)
    for i in eachindex(cellVolumeVal)
        v = cellVolumeVal[i]
        floor_v = floor(v)
        thick[i] = UInt16(min(Int(typemax(UInt16)),
                              Int(floor_v) + (mod(v, 1) > rand(rng) ? 1 : 0)))
    end
    idxAD = UInt16.(cellVolumeIdx) .* UInt16.(cellVolumeAD)
    idxBD = UInt16.(cellVolumeIdx) .* UInt16.(.!cellVolumeAD)
    valAD = thick .* UInt16.(cellVolumeAD)
    valBD = thick .* UInt16.(.!cellVolumeAD)

    _, dendnumAD = dilate_dendrite_paths_all(valAD, idxAD, neur_num; rng=rng)
    _, dendnumBD = dilate_dendrite_paths_all(valBD, idxBD, neur_num; rng=rng)

    for kk in 1:N_neur
        for x in gp_nuc[kk][1]
            dendnumAD[x] = 0
            dendnumBD[x] = 0
        end
        for x in gp_soma[kk]
            dendnumAD[x] = 0
            dendnumBD[x] = 0
        end
    end
    dendnumBD_merged = copy(dendnumBD)
    for i in eachindex(dendnumAD)
        if dendnumAD[i] > 0
            dendnumBD_merged[i] = dendnumAD[i]
        end
    end
    for i in eachindex(dendnumBD_merged)
        if dendnumBD_merged[i] > 0
            neur_num[i] = UInt16(dendnumBD_merged[i])
        end
    end
    for kk in 1:N_neur
        for x in gp_nuc[kk][1]
            neur_num[x] = 0
        end
        for x in gp_soma[kk]
            neur_num[x] = UInt16(kk)
        end
    end
    return neur_num, UInt16.(dendnumAD), gp_soma_aug
end

# ---------------------------------------------------------------------------
# Apical (through-volume) dendrites
# ---------------------------------------------------------------------------

"""
    grow_apical_dendrites!(vol_params, dend_params, neur_num, neur_num_AD_in,
                           gp_nuc, gp_soma; rng=Random.default_rng())
        -> (neur_num_out, neur_num_AD_out)

Add `N_den` through-volume apical dendrites rooted at random
unobstructed surface points. Single-stage Dijkstra; ports
`growApicalDendrites.m`.
"""
function grow_apical_dendrites!(vol_params::VolumeParams,
                                dend_params::DendriteParams,
                                neur_num::AbstractArray{<:Integer,3},
                                neur_num_AD_in::AbstractArray{<:Integer,3},
                                gp_nuc::AbstractVector,
                                gp_soma::AbstractVector;
                                rng::AbstractRNG=Random.default_rng())
    atParams       = collect(Float64, dend_params.atParams2)
    dweight        = dend_params.dweight
    bweight        = dend_params.bweight
    thicknessScale = dend_params.thicknessScale
    rallexp        = dend_params.rallexp
    vres           = vol_params.vres
    N_neur         = vol_params.N_neur
    N_den          = Int(round(vol_params.N_den))
    vol_sz         = vol_params.vol_sz
    fulldims = (Int(round(vol_sz[1] * vres)),
                Int(round(vol_sz[2] * vres)),
                Int(round(vol_sz[3] * vres)))
    atParams[2:4] .*= vres
    thicknessScale *= vres * vres

    cellVolume = Float32.(neur_num)
    for kk in 1:N_neur
        for x in gp_nuc[kk][1]
            cellVolume[x] = Float32(kk)
        end
    end
    cellVolumeIdx = zeros(Float32, fulldims)
    cellVolumeVal = zeros(Float32, fulldims)
    dendVar = 0.35

    roots = NTuple{3,Int}[]
    tries = 0
    max_tries = 1_000 * (N_den + 1)
    while length(roots) < N_den && tries < max_tries
        tries += 1
        x = max(1, ceil(Int, fulldims[1] * rand(rng)))
        y = max(1, ceil(Int, fulldims[2] * rand(rng)))
        if cellVolume[x, y, 1] == 0
            push!(roots, (x, y, fulldims[3]))
        end
    end
    lin_full = LinearIndices(fulldims)

    # Reusable Dijkstra buffers (see `_dijkstra!`); allocated once.
    dij_dist = fill(Inf32, fulldims)
    dij_path = zeros(Int, fulldims)
    dij_heap = _MinHeap()
    dij_touched = Int[]

    for j in 1:length(roots)
        rootL = roots[j]
        obstruction = copy(cellVolume)
        M = _build_dendrite_M(fulldims, obstruction, dweight, bweight; rng=rng)

        nat = Int(atParams[1])
        rootA = (clamp(Int(floor(rootL[1] + 2 * atParams[4] * (rand(rng) - 0.5))),
                       1, fulldims[1]),
                 clamp(Int(floor(rootL[2] + 2 * atParams[4] * (rand(rng) - 0.5))),
                       1, fulldims[2]),
                 fulldims[3])
        endsA = NTuple{3,Int}[]
        for _ in 1:nat
            flag = true
            distSC = 1.0
            numit = 0
            while flag && numit < 100
                θ = rand(rng) * 2π
                r = sqrt(rand(rng)) * atParams[2] * distSC
                x = clamp(Int(floor(r * cos(θ) + rootA[1])), 1, fulldims[1])
                y = clamp(Int(floor(r * sin(θ) + rootA[2])), 1, fulldims[2])
                z = 1
                if obstruction[x, y, z] == 0
                    push!(endsA, (x, y, z))
                    flag = false
                end
                distSC *= 1.01
                numit += 1
            end
            flag && push!(endsA, (rootA[1], rootA[2], 1))
        end

        # Endpoint sampling draws no Dijkstra output, so running it before
        # the search leaves the RNG stream untouched.
        targets = Int[lin_full[e...] for e in endsA]
        _dijkstra!(dij_dist, dij_path, dij_heap, dij_touched, M, rootL, targets)
        pathfrom = dij_path

        finepathsVal = zeros(Float32, fulldims)
        for e in endsA
            path = get_dendrite_path(pathfrom, e, rootL)
            isempty(path) && continue
            dendSz = max(0.0, 1 + dendVar * randn(rng))
            n = length(path)
            pw = ones(Float32, n)
            if n > 2
                for i in 2:(n - 1)
                    p, c, q = path[i - 1], path[i], path[i + 1]
                    d2 = abs(2c[1] - p[1] - q[1]) +
                         abs(2c[2] - p[2] - q[2]) +
                         abs(2c[3] - p[3] - q[3])
                    pw[i] = Float32(dendSz * (1 - (1 - 1 / sqrt(2)) * d2 / 2))
                end
                pw[1] = pw[2]; pw[end] = pw[end - 1]
                pw .= max.(pw, 0f0)
            else
                pw .= Float32(dendSz)
            end
            for (i, p) in enumerate(path)
                li = lin_full[p...]
                finepathsVal[li] += pw[i]
            end
        end
        for li in eachindex(finepathsVal)
            v = finepathsVal[li]
            v > 0 && (finepathsVal[li] = Float32(thicknessScale * atParams[5] *
                                                 v^(1 / rallexp)))
        end
        finepathsIdx = Float32(j + N_neur) .* Float32.(finepathsVal .> 0)
        cellVolume    .+= finepathsIdx
        cellVolumeIdx .+= finepathsIdx
        cellVolumeVal .+= finepathsVal
    end

    thick = UInt16.(min.(ceil.(cellVolumeVal), Float32(typemax(UInt16))))
    idx = UInt16.(cellVolumeIdx)
    _, dendnum = dilate_dendrite_paths_all(thick, idx, neur_num; rng=rng)
    cellVolumeAD_arr = UInt16.(neur_num_AD_in) .+ UInt16.(dendnum)
    neur_num_out = UInt16.(neur_num) .+ UInt16.(dendnum)
    for kk in 1:N_neur
        for x in gp_nuc[kk][1]
            neur_num_out[x] = 0
        end
    end
    for kk in 1:N_neur
        for x in gp_soma[kk]
            neur_num_out[x] = UInt16(kk)
        end
    end
    neur_num_AD_out = UInt16.(cellVolumeAD_arr)
    for kk in 1:N_neur
        for x in gp_nuc[kk][1]
            neur_num_AD_out[x] = 0
        end
        for x in gp_soma[kk]
            neur_num_AD_out[x] = 0
        end
    end
    for i in eachindex(neur_num_AD_out)
        if Int(neur_num_AD_out[i]) - Int(neur_num_out[i]) > 0
            neur_num_AD_out[i] = 0
        end
    end
    return neur_num_out, neur_num_AD_out
end

# ---------------------------------------------------------------------------
# Smooth soma-dendrite junction (simplified port)
# ---------------------------------------------------------------------------

"""
    smooth_cell_body(allpaths, cellBody, fdims) -> Vector{Int32}

Compute extra voxels that bridge dendrite roots and the soma body to
smooth the soma-dendrite junction. **Simplified port** of
`smoothCellBody.m`: rather than upstream's `cscvn` spline blend
between `connIdx` and `connRoots`, this adds a small radius-2 ball
around the first voxel where each path hits `cellBody`.
"""
function smooth_cell_body(allpaths::AbstractVector,
                          cellBody::AbstractVector{<:Integer},
                          fdims::NTuple{3,Int})
    cellBodySet = Set(Int32.(cellBody))
    H, W, D = fdims
    lin = LinearIndices((H, W, D))
    out = Set{Int32}()
    r = 2
    r2 = r * r
    for path in allpaths
        isempty(path) && continue
        hit = nothing
        for p in path
            li = Int32(lin[p...])
            if li in cellBodySet
                hit = p
                break
            end
        end
        hit === nothing && continue
        for dz in -r:r, dy in -r:r, dx in -r:r
            dx^2 + dy^2 + dz^2 <= r2 || continue
            ii = hit[1] + dx
            jj = hit[2] + dy
            kk = hit[3] + dz
            (1 <= ii <= H && 1 <= jj <= W && 1 <= kk <= D) || continue
            push!(out, Int32(lin[ii, jj, kk]))
        end
    end
    return collect(out)
end

# ---------------------------------------------------------------------------
# Per-cell fluorescence distribution
# ---------------------------------------------------------------------------

"""
    set_cell_fluorescence(vol_params, neur_params, dend_params,
                          neur_num, neur_soma, neur_num_AD, neur_locs, neur_vol;
                          rng=Random.default_rng())
        -> (gp_vals, neur_vol_out)

Populate per-voxel fluorescence values for each cell.

- **Soma** voxels: GP sample (via [`masked_3d_gp`](@ref)) normalised
  about a mean of 1.
- **Basal-dendrite** voxels: radial decay
  `wtSc[2]·exp(−r/(vres·wtSc[1])) + (1 − wtSc[2])` times `(1 − wtSc[3])`.
- **Apical-dendrite** voxels: uniform 1.
- **Background** processes (cells `N_neur+1 … N_neur+N_den`): uniform 1.

Returns `gp_vals` (per-cell `(loc, val, is_soma)` named tuples) plus
the updated `neur_vol`. Ports `setCellFluoresence.m`.
"""
function set_cell_fluorescence(vol_params::VolumeParams,
                               neur_params::NeuronParams,
                               dend_params::DendriteParams,
                               neur_num::AbstractArray{<:Integer,3},
                               neur_soma::AbstractArray{<:Integer,3},
                               neur_num_AD::AbstractArray{<:Integer,3},
                               neur_locs::AbstractMatrix{<:Real},
                               neur_vol::AbstractArray{<:Real,3};
                               rng::AbstractRNG=Random.default_rng())
    N_neur = vol_params.N_neur
    N_den  = Int(round(vol_params.N_den))
    vres   = vol_params.vres
    vol_sz = vol_params.vol_sz
    H = Int(round(vol_sz[1] * vres))
    W = Int(round(vol_sz[2] * vres))
    D = Int(round(vol_sz[3] * vres))
    cart = CartesianIndices((H, W, D))
    wtSc = collect(Float64, dend_params.weightScale)
    flSc = collect(Float64, neur_params.fluor_dist)

    numcomps = N_neur + N_den
    cell_voxels = [Int32[] for _ in 1:numcomps]
    for li in eachindex(neur_num)
        k = Int(neur_num[li])
        1 <= k <= numcomps && push!(cell_voxels[k], Int32(li))
    end

    gp_vals = Vector{NamedTuple{(:loc, :val, :is_soma),
                                Tuple{Vector{Int32}, Vector{Float32}, BitVector}}}(undef,
                                                                                   numcomps)
    neur_vol_out = Float32.(neur_vol)

    for kk in 1:N_neur
        loc = cell_voxels[kk]
        if isempty(loc)
            gp_vals[kk] = (loc=Int32[], val=Float32[], is_soma=BitVector(undef, 0))
            continue
        end
        is_soma_voxel = BitVector([neur_soma[Int(li)] == kk for li in loc])
        is_AD = BitVector([neur_num_AD[Int(li)] == kk for li in loc])

        gp_grid_size = max(1, Int(round(neur_params.avg_rad * 6 * vres)))
        gp = masked_3d_gp(gp_grid_size, flSc[1] * vres, flSc[2], 0.0; rng=rng)

        center_vox = (Int(floor(vres * neur_locs[kk, 1])),
                      Int(floor(vres * neur_locs[kk, 2])),
                      Int(floor(vres * neur_locs[kk, 3])))

        vals = zeros(Float32, length(loc))
        for (idx, li) in enumerate(loc)
            c = cart[Int(li)]
            r = sqrt((c[1] - vres * neur_locs[kk, 1])^2 +
                     (c[2] - vres * neur_locs[kk, 2])^2 +
                     (c[3] - vres * neur_locs[kk, 3])^2)
            base = Float32(wtSc[2] * exp(-r / (vres * wtSc[1])) +
                            (1 - wtSc[2])) * Float32(1 - wtSc[3])
            if is_soma_voxel[idx]
                gx = clamp(c[1] - center_vox[1] + cld(size(gp, 1), 2),
                           1, size(gp, 1))
                gy = clamp(c[2] - center_vox[2] + cld(size(gp, 2), 2),
                           1, size(gp, 2))
                gz = clamp(c[3] - center_vox[3] + cld(size(gp, 3), 2),
                           1, size(gp, 3))
                base = gp[gx, gy, gz]
            elseif is_AD[idx]
                base = 1f0
            end
            vals[idx] = base
        end
        si = findall(is_soma_voxel)
        if !isempty(si)
            m = sum(vals[si]) / length(si)
            scale = maximum(abs.(vals[si] .- m))
            if scale > 0
                vals[si] .= 0.5f0 .* (vals[si] .- m) ./ scale .+ 1f0
            else
                vals[si] .= 1f0
            end
            for i in si
                isnan(vals[i]) && (vals[i] = 1f0)
            end
        end
        gp_vals[kk] = (loc=loc, val=vals, is_soma=is_soma_voxel)
        for (idx, li) in enumerate(loc)
            neur_vol_out[Int(li)] = vals[idx]
        end
    end
    for kk in (N_neur + 1):numcomps
        loc = cell_voxels[kk]
        vals = ones(Float32, length(loc))
        gp_vals[kk] = (loc=loc, val=vals, is_soma=BitVector(falses(length(loc))))
        for li in loc
            neur_vol_out[Int(li)] = 1f0
        end
    end
    return gp_vals, neur_vol_out
end
