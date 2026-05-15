# Vasculature simulation.
#
# Ported from the upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - simulatebloodvessels.m → simulate_blood_vessels
#   - growMajorVessels.m     → grow_major_vessels!
#   - growCapillaries.m      → grow_capillaries!
#   - branchGrowNodes.m      → branch_grow_nodes!
#   - vessel_dijkstra.m      → vessel_dijkstra
#   - gennode.m / delnode.m  → VesselNode + gen_node / del_node!
#   - genconn.m              → VesselEdge + gen_conn
#   - nodesToConn.m          → nodes_to_conn
#   - connToVol.m            → conn_to_vol!
#   - pseudoRandSample2D.m   → pseudo_rand_sample_2d
#   - pseudoRandSample3D.m   → pseudo_rand_sample_3d

export VesselNode, VesselEdge,
       gen_node, gen_conn, del_node!, nodes_to_conn,
       pseudo_rand_sample_2d, pseudo_rand_sample_3d,
       vessel_dijkstra, branch_grow_nodes!,
       grow_major_vessels!, grow_capillaries!, conn_to_vol!,
       simulate_blood_vessels

# ---------------------------------------------------------------------------
# Node / edge primitives
# ---------------------------------------------------------------------------

"""
    VesselNode

A single branch point in the vasculature graph.

Mirrors upstream `gennode` struct fields. `root` uses two reserved values:
`0` means "no parent" (a source/root node) and `-1` means "deleted or
orphan" (upstream's empty-`[]` semantics; preserved so `delnode!` can blank
a slot without shifting indices). `type` selects the node role:
`:none, :edge, :surf, :vert, :sfvt, :capp`. `misc` is a generic per-node
payload (surface direction flags or leaf vessel-size).
"""
mutable struct VesselNode
    num::Int
    root::Int
    conn::Vector{Int}
    pos::Vector{Float64}
    type::Symbol
    misc::Vector{Float64}
end

"""
    gen_node()
    gen_node(num, root, conn, pos, type, misc=Float64[])

Construct a [`VesselNode`](@ref). The no-argument form returns a blank
sentinel used by [`del_node!`](@ref) to zero a slot without shifting
indices.
"""
gen_node() = VesselNode(0, -1, Int[], Float64[], :none, Float64[])

function gen_node(num::Integer, root, conn, pos, type::Symbol, misc=Float64[])
    rootv = isnothing(root) || (root isa AbstractVector && isempty(root)) ? -1 : Int(root)
    connv = isnothing(conn) || conn === () ? Int[] : collect(Int, conn)
    posv  = collect(Float64, pos)
    miscv = collect(Float64, misc)
    return VesselNode(Int(num), rootv, connv, posv, type, miscv)
end

"""
    VesselEdge

A single edge in the vasculature graph. Mirrors upstream `genconn` struct
fields. `weight = NaN` means "unset, fill on completion". `locs` is an
`N×3` matrix of voxel positions; empty (`0×3`) until [`conn_to_vol!`](@ref)
renders the edge. `misc` is a per-edge tag (`:none, :capp, :vtcp, ...`).
"""
mutable struct VesselEdge
    start::Int
    ends::Int
    weight::Float64
    locs::Matrix{Int}
    misc::Symbol
end

"""
    gen_conn()
    gen_conn(start, ends, weight=NaN, locs=zeros(Int,0,3), misc=:none)

Construct a [`VesselEdge`](@ref). No-argument form returns a blank sentinel.
"""
gen_conn() = VesselEdge(0, 0, NaN, zeros(Int, 0, 3), :none)

function gen_conn(start::Integer, ends::Integer, weight=NaN,
                  locs=zeros(Int, 0, 3), misc::Symbol=:none)
    w = isnothing(weight) ? NaN : Float64(weight)
    locm = locs === nothing || (locs isa AbstractVector && isempty(locs)) ?
        zeros(Int, 0, 3) : Matrix{Int}(locs)
    return VesselEdge(Int(start), Int(ends), w, locm, misc)
end

"""
    del_node!(nodes, num)

Mirrors upstream `delnode.m`. For every neighbour `i` in `nodes[num].conn`,
removes `num` from `nodes[i].conn`; if `nodes[i].root == num`, resets it to
the orphan sentinel `-1`. Finally replaces `nodes[num]` with a blank
[`VesselNode`](@ref), preserving array indices.
"""
function del_node!(nodes::AbstractVector{VesselNode}, num::Integer)
    for i in nodes[num].conn
        filter!(c -> c != num, nodes[i].conn)
        if nodes[i].root == num
            nodes[i].root = -1
        end
    end
    nodes[num] = gen_node()
    return nodes
end

"""
    nodes_to_conn(nodes) -> Vector{VesselEdge}

Construct edges from the node tree by walking from each leaf back through
`root` pointers. Weights are accumulated by `sqrt(w² + leaf_misc²)` along
the path (matching upstream `nodesToConn.m`).
"""
function nodes_to_conn(nodes::AbstractVector{VesselNode})
    N = length(nodes)
    # weight matrix, sparse via dict (only edges that exist)
    W = Dict{Tuple{Int,Int}, Float64}()
    leaf_indices = [i for i in 1:N if length(nodes[i].conn) == 1]
    for li in leaf_indices
        leaf_size = isempty(nodes[li].misc) ? 0.0 : nodes[li].misc[1]
        cur = li
        while nodes[cur].root > 0
            r = nodes[cur].root
            key = (nodes[cur].num, r)
            prev = get(W, key, 0.0)
            W[key] = sqrt(prev^2 + leaf_size^2)
            cur = r
        end
    end
    conn = VesselEdge[]
    for ((s, e), w) in W
        push!(conn, gen_conn(s, e, w))
    end
    return conn
end

# ---------------------------------------------------------------------------
# Pseudo-random pseudo-uniform spatial sampling
# ---------------------------------------------------------------------------

function _correct_offset_round(c::AbstractVector, i::AbstractVector)
    # Equivalent of upstream's correctOffsetRoundError: keep window spans equal.
    dc = c[2] - c[1]
    di = i[2] - i[1]
    if dc > di
        return (c[1] + (dc - di), c[2]), (i[1], i[2])
    elseif dc < di
        return (c[1], c[2]), (i[1] + (di - dc), i[2])
    else
        return (c[1], c[2]), (i[1], i[2])
    end
end

"""
    pseudo_rand_sample_2d(sz, nsamps; width=2.0, weight=1.0, pdf=ones(sz), maxit=1000, rng=Random.default_rng())

Rejection-sample `nsamps` positions in an `sz` grid with a Gaussian
exclusionary kernel (width `width`, depth `weight`) around each accepted
sample. Returns `pos::Matrix{Int}` of size `nsamps × 2` and the updated
`pdf`. Ports `pseudoRandSample2D.m`.
"""
function pseudo_rand_sample_2d(sz::NTuple{2,<:Integer}, nsamps::Integer;
                               width::Real=2.0, weight::Real=1.0,
                               pdf::AbstractMatrix{<:Real}=ones(Float32, sz),
                               maxit::Integer=1000,
                               rng::AbstractRNG=Random.default_rng())
    pdf = copy(Matrix{Float32}(pdf))
    w2 = Int(ceil(2 * width))
    gpdf = Matrix{Float32}(undef, 2 * w2 + 1, 2 * w2 + 1)
    for ix in -w2:w2, iy in -w2:w2
        gpdf[ix + w2 + 1, iy + w2 + 1] = -Float32(weight) * exp(-(ix^2 + iy^2) / width^2)
    end
    pos = zeros(Int, nsamps, 2)
    i = 1
    numit = 0
    while i <= nsamps && numit < maxit
        numit += 1
        rx = max(1, ceil(Int, rand(rng) * sz[1]))
        ry = max(1, ceil(Int, rand(rng) * sz[2]))
        if pdf[rx, ry] - rand(rng) > 0
            xc = (max(0, w2 + 1 - rx) + 1, min(0, sz[1] - rx - w2) + 2 * w2 + 1)
            yc = (max(0, w2 + 1 - ry) + 1, min(0, sz[2] - ry - w2) + 2 * w2 + 1)
            xi = (max(1, rx - w2), min(sz[1], rx + w2))
            yi = (max(1, ry - w2), min(sz[2], ry + w2))
            xc, xi = _correct_offset_round(collect(xc), collect(xi))
            yc, yi = _correct_offset_round(collect(yc), collect(yi))
            @views pdf[xi[1]:xi[2], yi[1]:yi[2]] .+= gpdf[xc[1]:xc[2], yc[1]:yc[2]]
            pos[i, 1] = rx
            pos[i, 2] = ry
            i += 1
            numit = 0
        end
    end
    if i <= nsamps
        pos = pos[1:i - 1, :]
    end
    return pos, pdf
end

"""
    pseudo_rand_sample_3d(sz, nsamps; width=2.0, weight=1.0, pdf=ones(sz), maxit=10_000, rng=Random.default_rng())

3-D analogue of [`pseudo_rand_sample_2d`](@ref). Returns `pos::Matrix{Int}`
of size `nsamps × 3` and the updated `pdf`.
"""
function pseudo_rand_sample_3d(sz::NTuple{3,<:Integer}, nsamps::Integer;
                               width::Real=2.0, weight::Real=1.0,
                               pdf::AbstractArray{<:Real,3}=ones(Float32, sz),
                               maxit::Integer=10_000,
                               rng::AbstractRNG=Random.default_rng())
    pdf = copy(Array{Float32,3}(pdf))
    w2 = Int(ceil(2 * width))
    gpdf = Array{Float32,3}(undef, 2 * w2 + 1, 2 * w2 + 1, 2 * w2 + 1)
    for ix in -w2:w2, iy in -w2:w2, iz in -w2:w2
        gpdf[ix + w2 + 1, iy + w2 + 1, iz + w2 + 1] =
            -Float32(weight) * exp(-(ix^2 + iy^2 + iz^2) / width^2)
    end
    pos = zeros(Int, nsamps, 3)
    i = 1
    numit = 0
    while i <= nsamps && numit < maxit
        numit += 1
        rx = max(1, ceil(Int, rand(rng) * sz[1]))
        ry = max(1, ceil(Int, rand(rng) * sz[2]))
        rz = max(1, ceil(Int, rand(rng) * sz[3]))
        if pdf[rx, ry, rz] - rand(rng) > 0
            xc = (max(0, w2 + 1 - rx) + 1, min(0, sz[1] - rx - w2) + 2 * w2 + 1)
            yc = (max(0, w2 + 1 - ry) + 1, min(0, sz[2] - ry - w2) + 2 * w2 + 1)
            zc = (max(0, w2 + 1 - rz) + 1, min(0, sz[3] - rz - w2) + 2 * w2 + 1)
            xi = (max(1, rx - w2), min(sz[1], rx + w2))
            yi = (max(1, ry - w2), min(sz[2], ry + w2))
            zi = (max(1, rz - w2), min(sz[3], rz + w2))
            xc, xi = _correct_offset_round(collect(xc), collect(xi))
            yc, yi = _correct_offset_round(collect(yc), collect(yi))
            zc, zi = _correct_offset_round(collect(zc), collect(zi))
            @views pdf[xi[1]:xi[2], yi[1]:yi[2], zi[1]:zi[2]] .+=
                gpdf[xc[1]:xc[2], yc[1]:yc[2], zc[1]:zc[2]]
            pos[i, 1] = rx
            pos[i, 2] = ry
            pos[i, 3] = rz
            i += 1
            numit = 0
        end
    end
    if i <= nsamps
        pos = pos[1:i - 1, :]
    end
    return pos, pdf
end

# ---------------------------------------------------------------------------
# Binary morphology (small, hand-rolled — kept off ImageMorphology for now)
# ---------------------------------------------------------------------------

function _disk_strel_2d(r::Real)
    rr = Int(ceil(r))
    mask = falses(2 * rr + 1, 2 * rr + 1)
    r2 = r^2
    for i in -rr:rr, j in -rr:rr
        mask[i + rr + 1, j + rr + 1] = i^2 + j^2 <= r2
    end
    return mask
end

"""
    dilate2d_disk!(mask::AbstractMatrix{Bool}, r::Real)

Binary 2-D dilation by a disk of radius `r` (MATLAB
`imdilate(mask, strel('disk', r))`). Operates in-place; returns `mask`.
"""
function dilate2d_disk!(mask::AbstractMatrix{Bool}, r::Real)
    rr = Int(ceil(r))
    rr == 0 && return mask
    se = _disk_strel_2d(r)
    src = copy(mask)
    fill!(mask, false)
    H, W = size(mask)
    for j in 1:W, i in 1:H
        src[i, j] || continue
        for dj in -rr:rr, di in -rr:rr
            se[di + rr + 1, dj + rr + 1] || continue
            ii = i + di
            jj = j + dj
            (1 <= ii <= H && 1 <= jj <= W) || continue
            mask[ii, jj] = true
        end
    end
    return mask
end

"""
    paint_ball3d!(vol, p::AbstractVector{<:Integer}, r::Real, val=true)

Set every voxel within Euclidean radius `r` of `p` (1-based integer
coordinates) to `val`, clipping at array bounds. Equivalent to
`imdilate(point_mask, strel(ball(r)))` for a single point.
"""
function paint_ball3d!(vol::AbstractArray{T,3}, p::AbstractVector{<:Integer},
                       r::Real, val::T=true) where {T}
    rr = Int(ceil(r))
    rr == 0 && (vol[p[1], p[2], p[3]] = val; return vol)
    r2 = r^2
    H, W, D = size(vol)
    for dk in -rr:rr, dj in -rr:rr, di in -rr:rr
        di^2 + dj^2 + dk^2 <= r2 || continue
        ii = p[1] + di
        jj = p[2] + dj
        kk = p[3] + dk
        (1 <= ii <= H && 1 <= jj <= W && 1 <= kk <= D) || continue
        vol[ii, jj, kk] = val
    end
    return vol
end

# ---------------------------------------------------------------------------
# Dijkstra (dense distance matrix; Inf marks forbidden)
# ---------------------------------------------------------------------------

"""
    vessel_dijkstra(distMat, proot) -> (distance, pathfrom)

Direct port of upstream `vessel_dijkstra.m`. Operates on a dense
distance matrix where `Inf` marks a forbidden edge. `distance[i]` is the
minimum path length from `proot` to `i` (`Inf` if unreachable);
`pathfrom[i]` is the predecessor of `i` on that path (`0` for the root
and for unreachable nodes).
"""
function vessel_dijkstra(distMat::AbstractMatrix{<:Real}, proot::Integer)
    N = size(distMat, 1)
    distance = fill(Inf, N)
    pathfrom = zeros(Int, N)
    N == 0 && return distance, pathfrom
    tovisit  = trues(N)
    unvisited = zeros(Float64, N)
    unvisited[proot] = 1.0      # sentinel >0 so it gets popped first
    distance[proot] = 0.0
    cn = Int(proot)

    while any(!iszero, unvisited)
        tovisit[cn] = false
        # pick the smallest *nonzero* entry of unvisited
        best = Inf
        cn_next = 0
        @inbounds for k in 1:N
            u = unvisited[k]
            if u != 0 && u < best
                best = u
                cn_next = k
            end
        end
        cn = cn_next
        unvisited[cn] = 0.0
        d_cn = distance[cn]
        @inbounds for nn in 1:N
            tovisit[nn] || continue
            ndist = d_cn + distMat[cn, nn]
            if ndist < distance[nn]
                unvisited[nn] = max(eps(), ndist)
                distance[nn]  = ndist
                pathfrom[nn]  = cn
            end
        end
    end
    return distance, pathfrom
end

# ---------------------------------------------------------------------------
# Surface vessel growth
# ---------------------------------------------------------------------------

function _veclinspace(a::AbstractVector, b::AbstractVector, n::Integer)
    n = max(2, Int(n))
    out = zeros(Int, 2, n)
    for k in 1:n
        t = (k - 1) / (n - 1)
        out[1, k] = Int(round(a[1] + (b[1] - a[1]) * t))
        out[2, k] = Int(round(a[2] + (b[2] - a[2]) * t))
    end
    return out
end

"""
    branch_grow_nodes!(nodes, neur_ves, params, idx, dir; rng=Random.default_rng())

Grow surface vasculature outward from node `idx` along direction `dir`
(radians). Avoids overlapping previously placed paths. Mutates `nodes`
(appending new entries) and `neur_ves` (a 2-D occupancy mask). Ports
`branchGrowNodes.m`.
"""
function branch_grow_nodes!(nodes::Vector{VesselNode},
                            neur_ves::AbstractMatrix{Bool},
                            params::VasculatureNodeParams,
                            idx::Integer,
                            dir::Real;
                            rng::AbstractRNG=Random.default_rng())
    borderflag  = true
    overlapflag = false
    numIt       = 0
    prevPos     = nodes[idx].pos[1:2]
    prevNum     = nodes[idx].num
    testIdxs2   = Tuple{Int,Int}[]
    nv_size     = size(neur_ves)
    branchP     = 0.0
    vesrad      = max(1, Int(ceil(params.vesrad)))
    neur_ves2   = copy(neur_ves)
    dilate2d_disk!(neur_ves2, vesrad)

    while numIt < params.maxit && borderflag
        if rand(rng) < branchP
            branchP = 0.0
            dirB = dir - abs((0.5 + rand(rng)) * params.dirvar)
            dirA = dir + abs((0.5 + rand(rng)) * params.dirvar)
            dir = dirA
            branch_grow_nodes!(nodes, neur_ves, params, prevNum, dirB; rng)
            overlapflag = true
        else
            branchP += params.branchp
        end

        dirVect = (cos(dir), sin(dir))
        vesDist = max(params.varsc * randn(rng) + params.lensc, params.mindist)
        nodePos = [dirVect[1] * vesDist + prevPos[1] + params.varpos * randn(rng),
                   dirVect[2] * vesDist + prevPos[2] + params.varpos * randn(rng)]
        if nodePos[1] < 1 || nodePos[2] < 1 ||
           nodePos[1] > nv_size[1] || nodePos[2] > nv_size[2]
            nodePos[1] = clamp(nodePos[1], 1.0, Float64(nv_size[1]))
            nodePos[2] = clamp(nodePos[2], 1.0, Float64(nv_size[2]))
            borderflag = false
        end

        ts = _veclinspace(nodePos, prevPos, Int(ceil(vesDist)))
        testIdxs = Tuple{Int,Int}[]
        for k in 1:size(ts, 2)
            for (di, dj) in ((0,0), (0,1), (1,0))
                ii = clamp(ts[1, k] + di, 1, nv_size[1])
                jj = clamp(ts[2, k] + dj, 1, nv_size[2])
                push!(testIdxs, (ii, jj))
            end
        end

        hit = false
        if !overlapflag
            for (ii, jj) in testIdxs
                if neur_ves2[ii, jj]
                    hit = true
                    break
                end
            end
        end

        if !hit
            nodeNum = length(nodes) + 1
            push!(nodes, gen_node(nodeNum, prevNum, [prevNum], nodePos, :surf))
            prevPos = nodePos
            prevNum = nodeNum
            numIt += 1
            for (ii, jj) in testIdxs2
                neur_ves[ii, jj] = true
            end
            testIdxs2 = testIdxs
        else
            borderflag = false
            testIdxs   = Tuple{Int,Int}[]
            testIdxs2  = Tuple{Int,Int}[]
        end
        overlapflag = false
    end
    for (ii, jj) in testIdxs2
        neur_ves[ii, jj] = true
    end
    return nodes, neur_ves
end

# ---------------------------------------------------------------------------
# Major-vessel growth (surface + diving)
# ---------------------------------------------------------------------------

function _pos2dists(pos::AbstractMatrix{<:Real})
    N = size(pos, 1)
    D = zeros(Float64, N, N)
    for j in 1:N, i in 1:N
        if i != j
            D[i, j] = sqrt((pos[i, 1] - pos[j, 1])^2 +
                           (pos[i, 2] - pos[j, 2])^2 +
                           (pos[i, 3] - pos[j, 3])^2)
        end
    end
    return D
end

"""
    grow_major_vessels!(np::VasculatureNodeParams, vp::NamedTuple, nv::NamedTuple; rng=Random.default_rng())
        -> (nodes, nv_updated)

Generate source/edge/surface/diving nodes for the major-vessel graph.
`vp` carries the scaled vasculature parameters (`depth_surf`, `vesSize`,
`mindists`, `sepweight`, `distWeightScale`, `randWeightScale`); `nv` carries
volume-size and target counts. Ports `growMajorVessels.m`.
"""
function grow_major_vessels!(np::VasculatureNodeParams, vp::NamedTuple,
                             nv::NamedTuple;
                             rng::AbstractRNG=Random.default_rng())
    sz = nv.size
    nodes = VesselNode[]
    nsource = nv.nsource
    nsurf   = nv.nsurf

    # Seed source nodes on the four edges of the surface plane (z = depth_surf).
    for i in 1:nsource
        b1 = rand(rng) >= (sz[2] / (sz[1] + sz[2]))
        b2 = rand(rng) >= 0.5
        if      b1 &&  b2
            TMPPOS = [max(1, ceil(Int, sz[1] * rand(rng))), 1, vp.depth_surf]
        elseif  b1 && !b2
            TMPPOS = [max(1, ceil(Int, sz[1] * rand(rng))), sz[2], vp.depth_surf]
        elseif !b1 &&  b2
            TMPPOS = [1, max(1, ceil(Int, sz[2] * rand(rng))), vp.depth_surf]
        else
            TMPPOS = [sz[1], max(1, ceil(Int, sz[2] * rand(rng))), vp.depth_surf]
        end
        push!(nodes, gen_node(i, 0, Int[], TMPPOS, :edge, [b1 ? 1.0 : 0.0, b2 ? 1.0 : 0.0]))
    end

    neur_surf = falses(sz[1], sz[2])
    for i in 1:nsource
        b1 = nodes[i].misc[1] == 1.0
        b2 = nodes[i].misc[2] == 1.0
        randDir = if      b1 &&  b2; 0.5π + randn(rng) * np.dirvar
                  elseif  b1 && !b2; 1.5π + randn(rng) * np.dirvar
                  elseif !b1 &&  b2; 0.0π + randn(rng) * np.dirvar
                  else                 1.0π + randn(rng) * np.dirvar
                  end
        branch_grow_nodes!(nodes, neur_surf, np, i, randDir; rng)
    end
    nlinks = length(nodes) - nsource

    # Fill in z=depth_surf for newly grown 2-D surface nodes.
    for i in (nsource + 1):(nlinks + nsource)
        nodes[i].pos = [round(nodes[i].pos[1]),
                        round(nodes[i].pos[2]),
                        Float64(vp.depth_surf)]
    end

    # Dilate the surface occupancy and sample additional surf vertices.
    dilate2d_disk!(neur_surf, 2 * max(1, Int(ceil(np.vesrad))))
    pdf0 = Float32.(1 .- neur_surf)
    surfpos_raw, _ = pseudo_rand_sample_2d(sz[1:2], nsurf;
                                           width=vp.mindists[1],
                                           weight=vp.sepweight,
                                           pdf=pdf0,
                                           rng=rng)
    actual_nsurf = size(surfpos_raw, 1)
    surfpos_new = hcat(surfpos_raw, fill(Int(vp.depth_surf), actual_nsurf))

    # Stack existing edge/surface node positions on top of the new samples.
    npre = nlinks + nsource
    existing = zeros(Float64, npre, 3)
    for i in 1:npre
        existing[i, :] = nodes[i].pos
    end
    surfpos = vcat(existing, Float64.(surfpos_new))

    # Distance matrix with structural Inf/0 blocks per upstream.
    surfmat = _pos2dists(surfpos)
    surfmat[1:npre, 1:npre] .= Inf
    surfmat[1:nsource, 1:nsource] .= 0.0
    for i in 1:npre
        r = nodes[i].root
        if r > 0
            d = sqrt(sum((nodes[i].pos .- nodes[r].pos).^2))
            surfmat[r, i] = d
            surfmat[i, r] = d
        end
    end
    # Forbid edges *from* the new surf points back to existing nodes (only
    # forward direction is allowed).
    surfmat[(npre + 1):end, 1:npre] .= Inf

    surfpath = zeros(Int, size(surfmat, 1))
    if size(surfmat, 1) > 0 && nsource > 0
        TMPsurfmat = (surfmat .^ vp.distWeightScale) .*
                     (1 .+ vp.randWeightScale .* randn(rng, size(surfmat)))
        _, surfpath = vessel_dijkstra(TMPsurfmat, 1)
        surfpath[1:nsource] .= 1:nsource
    end
    for i in 1:min(actual_nsurf + nsource, length(surfpath))
        if surfpath[i] == 1
            best = 1
            best_d = Inf
            for k in 1:nsource
                if surfmat[i, k] < best_d
                    best_d = surfmat[i, k]
                    best = k
                end
            end
            surfpath[i] = best
        end
    end

    # Insert the new surf nodes.
    for j in 1:actual_nsurf
        i = npre + j
        sp = i <= length(surfpath) ? surfpath[i] : 0
        is_edge = surfpos[i, 1] == 1 || surfpos[i, 1] == sz[1] ||
                  surfpos[i, 2] == 1 || surfpos[i, 2] == sz[2]
        root = sp > 0 ? sp : 0
        connlist = sp > 0 ? [sp] : Int[]
        push!(nodes,
              gen_node(i, root, connlist, surfpos[i, :], is_edge ? :edge : :surf))
    end
    nnodes = npre + actual_nsurf

    # Backfill conn lists from root pointers.
    for i in 1:nnodes
        r = nodes[i].root
        if r > 0
            if !(i in nodes[r].conn)
                push!(nodes[r].conn, i)
                sort!(nodes[r].conn)
            end
        end
    end

    # Prune lone surf nodes and elect some surf nodes as diving (sfvt) nodes.
    se_rad = max(1, Int(round(vp.mindists[1] * 2)))
    neur_vert = falses(sz[1], sz[2])
    for i in 1:nnodes
        if nodes[i].type === :surf && length(nodes[i].conn) == 1
            x, y = Int(round(nodes[i].pos[1])), Int(round(nodes[i].pos[2]))
            if neur_vert[x, y] == false
                nodes[i].type = :sfvt
                tmp = falses(sz[1], sz[2])
                tmp[x, y] = true
                dilate2d_disk!(tmp, se_rad)
                neur_vert .|= tmp
            else
                del_node!(nodes, i)
            end
        end
    end

    surfidx = [i for i in 1:nnodes if nodes[i].type === :surf]
    while count(n -> n.type === :sfvt, nodes) < nv.nvert
        candidates = [i for i in surfidx
                      if !neur_vert[Int(round(nodes[i].pos[1])),
                                    Int(round(nodes[i].pos[2]))]]
        isempty(candidates) && break
        TMPIDX = candidates[max(1, ceil(Int, rand(rng) * length(candidates)))]
        nodes[TMPIDX].type = :sfvt
        tmp = falses(sz[1], sz[2])
        x, y = Int(round(nodes[TMPIDX].pos[1])), Int(round(nodes[TMPIDX].pos[2]))
        tmp[x, y] = true
        dilate2d_disk!(tmp, se_rad)
        neur_vert .|= tmp
    end

    # Grow diving vessels down to the bottom of the volume.
    vertidx = [i for i in 1:nnodes if nodes[i].type === :sfvt]
    TMPIDX = nnodes
    last_top = TMPIDX
    for vi in vertidx
        curr_node = vi
        while nodes[curr_node].pos[3] < sz[3]
            TMPIDX += 1
            zinc = max(np.varsc * randn(rng) + np.lensc, np.mindist)
            np_pos = [nodes[curr_node].pos[1] + ceil(randn(rng) * np.varpos),
                      nodes[curr_node].pos[2] + ceil(randn(rng) * np.varpos),
                      nodes[curr_node].pos[3] + ceil(zinc)]
            np_pos[1] = clamp(np_pos[1], 1.0, Float64(sz[1]))
            np_pos[2] = clamp(np_pos[2], 1.0, Float64(sz[2]))
            np_pos[3] = clamp(np_pos[3], 1.0, Float64(sz[3]))
            push!(nodes, gen_node(TMPIDX, curr_node, [curr_node], np_pos, :vert))
            r = nodes[TMPIDX].root
            if !(TMPIDX in nodes[r].conn)
                push!(nodes[r].conn, TMPIDX)
                sort!(nodes[r].conn)
            end
            curr_node = TMPIDX
        end
        last_top = TMPIDX
    end
    nvert    = count(n -> n.type === :sfvt, nodes)
    nvertconn = last_top - nnodes
    nnodes = last_top

    # End nodes: gamma-distributed vessel sizes.
    ends_idx = [i for i in 1:nnodes if length(nodes[i].conn) == 1]
    γshape = 3.0
    γscale = (vp.vesSize[2] - vp.vesSize[3]) / 3
    for ei in ends_idx
        s = γscale > 0 ? rand(rng, Gamma(γshape, γscale)) : 0.0
        nodes[ei].misc = [vp.vesSize[3] + s]
    end

    nv2 = merge(nv, (nlinks = nlinks, nnodes = nnodes,
                     nvert = nvert, nvertconn = nvertconn))
    return nodes, nv2
end

# ---------------------------------------------------------------------------
# Capillary growth
# ---------------------------------------------------------------------------

"""
    grow_capillaries!(nodes, conn, neur_ves, nv::NamedTuple, vp::NamedTuple, vres; rng=Random.default_rng())
        -> (nodes, conn, nv_updated)

Generate capillary nodes and capillary-to-capillary / capillary-to-diving
connections. Ports `growCapillaries.m`.
"""
function grow_capillaries!(nodes::Vector{VesselNode},
                           conn::Vector{VesselEdge},
                           neur_ves::AbstractArray{Bool,3},
                           nv::NamedTuple, vp::NamedTuple,
                           vres::Real;
                           rng::AbstractRNG=Random.default_rng())
    szum = nv.szum
    szum = (Int(szum[1]), Int(szum[2]), Int(szum[3]))
    sz = nv.size

    # Build a smoothed exclusion volume from existing connection locations.
    dilrad = max(1, Int(ceil(vp.mindists[1] / vres)))
    se = Array{Float32,3}(undef, 2 * dilrad + 1, 2 * dilrad + 1, 2 * dilrad + 1)
    for ix in -dilrad:dilrad, iy in -dilrad:dilrad, iz in -dilrad:dilrad
        se[ix + dilrad + 1, iy + dilrad + 1, iz + dilrad + 1] =
            exp(-2 * (ix^2 + iy^2 + iz^2) / dilrad^2)
    end
    TMPvol = zeros(Float32, szum)
    step = max(1, Int(floor(dilrad / 3)))
    for c in conn
        isempty(c.locs) && continue
        Lrows = c.locs
        for ridx in 1:step:size(Lrows, 1)
            p = Int.(ceil.(Lrows[ridx, :] ./ vres))
            for ax in 1:3
                p[ax] = clamp(p[ax], 1, szum[ax])
            end
            for di in -dilrad:dilrad, dj in -dilrad:dilrad, dk in -dilrad:dilrad
                ii = p[1] + di
                jj = p[2] + dj
                kk = p[3] + dk
                (1 <= ii <= szum[1] && 1 <= jj <= szum[2] && 1 <= kk <= szum[3]) || continue
                v = se[di + dilrad + 1, dj + dilrad + 1, dk + dilrad + 1]
                if v > TMPvol[ii, jj, kk]
                    TMPvol[ii, jj, kk] = v
                end
            end
        end
    end

    capppos_um, _ = pseudo_rand_sample_3d(szum, nv.ncapp;
                                          width=vp.mindists[3] / vres,
                                          weight=vp.sepweight,
                                          pdf=1.0f0 .- TMPvol,
                                          rng=rng)
    actual_ncapp = size(capppos_um, 1)
    capppos = zeros(Int, actual_ncapp, 3)
    @inbounds for i in 1:actual_ncapp, k in 1:3
        capppos[i, k] = Int(capppos_um[i, k]) * Int(vres) + 1 -
                        max(1, ceil(Int, rand(rng) * max(1, Int(vres))))
        capppos[i, k] = clamp(capppos[i, k], 1, sz[k])
    end

    # Decide how many capp-to-vert connections per diving vessel.
    nvert_conn = Int[max(1, ceil(Int, rand(rng) * max(1, ceil(Int, szum[3] / vp.vesFreq[3]))))
                     for _ in 1:nv.nvert]
    nvert_sum = isempty(nvert_conn) ? 0 : sum(nvert_conn)
    nodeIdx = nv.nnodes
    connIdx = length(conn)
    vesSize = length(vp.vesSize) >= 4 ? vp.vesSize : vcat(vp.vesSize, 0.0)

    vertidxs = [i for i in 1:nv.nnodes if nodes[i].type === :sfvt]
    capp_active = trues(actual_ncapp)
    for (i, vesidx_start) in enumerate(vertidxs)
        i > length(nvert_conn) && break
        # Build the chain of diving-vessel descendants reachable from vesidx_start.
        chain = Int[vesidx_start]
        flag = true
        while flag
            last = chain[end]
            kids = [c for c in nodes[last].conn if c != last &&
                    (1 <= c <= length(nodes)) && nodes[c].type === :vert]
            kids = setdiff(kids, chain)
            if isempty(kids)
                flag = false
            else
                push!(chain, first(kids))
            end
        end
        vesidx = chain[2:end]
        isempty(vesidx) && continue
        for _ in 1:nvert_conn[i]
            TMPidx = vesidx[max(1, ceil(Int, rand(rng) * length(vesidx)))]
            nodeIdx += 1
            connIdx += 1
            # Find nearest active capp position.
            tgt = nodes[TMPidx].pos
            best = 0
            best_d = Inf
            @inbounds for k in 1:actual_ncapp
                capp_active[k] || continue
                d = (capppos[k, 1] - tgt[1])^2 +
                    (capppos[k, 2] - tgt[2])^2 +
                    (capppos[k, 3] - tgt[3])^2
                if d < best_d
                    best_d = d
                    best = k
                end
            end
            best == 0 && break
            cap_xyz = Float64.(capppos[best, :])
            push!(nodes, gen_node(nodeIdx, TMPidx, [TMPidx], cap_xyz, :capp))
            if !(nodeIdx in nodes[TMPidx].conn)
                push!(nodes[TMPidx].conn, nodeIdx)
                sort!(nodes[TMPidx].conn)
            end
            w = max(1.0, vesSize[3] + vesSize[4] * randn(rng))
            push!(conn, gen_conn(nodeIdx, TMPidx, w, zeros(Int, 0, 3), :vtcp))
            capp_active[best] = false
        end
    end
    nv = merge(nv, (nnodes = nodeIdx, nconn = connIdx, nvert_sum = nvert_sum))

    # Remaining capillaries become orphan capp nodes (root = -1).
    for k in 1:actual_ncapp
        capp_active[k] || continue
        nodeIdx += 1
        push!(nodes, gen_node(nodeIdx, -1, Int[], Float64.(capppos[k, :]), :capp))
    end
    nv = merge(nv, (nnodes = nodeIdx,))

    # Capillary connectivity matrix.
    vert_capp_idxs = [i for i in 1:nv.nnodes
                      if nodes[i].type === :capp && nodes[i].root >= 0]
    orphan_capp_idxs = [i for i in 1:nv.nnodes
                        if nodes[i].type === :capp && nodes[i].root < 0]
    connidxs = vcat(vert_capp_idxs, orphan_capp_idxs)
    ncapp = length(connidxs)
    capppos2 = zeros(Float64, ncapp, 3)
    for (k, ix) in enumerate(connidxs)
        capppos2[k, :] = nodes[ix].pos
    end

    cappmat = _pos2dists(capppos2)
    for k in 1:ncapp; cappmat[k, k] = Inf; end
    cappmat[1:nvert_sum, 1:nvert_sum] .= Inf
    cappconnmat = falses(ncapp, ncapp)
    if ncapp > 1
        @inbounds for k in 1:ncapp
            best = 1
            best_d = Inf
            for j in 1:ncapp
                if cappmat[k, j] < best_d
                    best_d = cappmat[k, j]
                    best = j
                end
            end
            cappconnmat[k, best] = true
            cappconnmat[best, k] = true
            cappmat[k, best] = Inf
            cappmat[best, k] = Inf
        end
    end
    @inbounds for k in 1:ncapp, j in 1:ncapp
        if cappmat[k, j] > vp.maxcappdist
            cappmat[k, j] = Inf
        end
    end
    for k in 1:ncapp
        if sum(cappconnmat[k, :]) >= 3
            cappmat[k, :] .= Inf
            cappmat[:, k] .= Inf
        end
        for j in 1:ncapp
            if cappconnmat[k, j]
                cappmat[k, j] = Inf
                cappmat[j, k] = Inf
            end
        end
    end

    # Vessel-mask blocking: forbid pairs whose line crosses a vessel.
    H, W, D = size(neur_ves)
    for i in (nvert_sum + 1):ncapp
        for j in (i + 1):ncapp
            cappmat[i, j] < Inf || continue
            nseg = max(2, 2 * Int(ceil(cappmat[i, j])))
            blocked = false
            @inbounds for t in 0:(nseg - 1)
                α = t / (nseg - 1)
                xi = Int(ceil((1 - α) * capppos2[i, 1] + α * capppos2[j, 1]))
                yi = Int(ceil((1 - α) * capppos2[i, 2] + α * capppos2[j, 2]))
                zi = Int(ceil((1 - α) * capppos2[i, 3] + α * capppos2[j, 3]))
                xi = clamp(xi, 1, H)
                yi = clamp(yi, 1, W)
                zi = clamp(zi, 1, D)
                if neur_ves[xi, yi, zi]
                    blocked = true
                    break
                end
            end
            if blocked
                cappmat[i, j] = Inf
                cappmat[j, i] = Inf
            end
        end
    end

    # Iteratively connect under-degree capp nodes.
    safety = 100 * ncapp + 1
    while safety > 0
        safety -= 1
        cappsum = vec(sum(cappconnmat, dims=2))
        if (ncapp - nvert_sum > 0) && all(c -> c > 1, @view(cappsum[(nvert_sum + 1):end]))
            break
        end
        idxs = findall(==(1), cappsum)
        isempty(idxs) && break
        if all(>=(Inf), @view cappmat[:, idxs])
            break
        end
        rndidx = idxs[max(1, ceil(Int, rand(rng) * length(idxs)))]
        if minimum(@view cappmat[rndidx, :]) < Inf
            inv_d = 1.0 ./ (cappmat[rndidx, :] .^ vp.distsc)
            S = sum(inv_d)
            S > 0 || break
            cdf = cumsum(inv_d) ./ S
            u = rand(rng)
            lnkidx = findfirst(>=(u), cdf)
            lnkidx === nothing && break
            cappconnmat[rndidx, lnkidx] = true
            cappconnmat[lnkidx, rndidx] = true
            for j in 1:ncapp
                if cappconnmat[rndidx, j]
                    cappmat[j, lnkidx] = Inf
                    cappmat[lnkidx, j] = Inf
                end
                if cappconnmat[lnkidx, j]
                    cappmat[j, rndidx] = Inf
                    cappmat[rndidx, j] = Inf
                end
            end
            cappmat[rndidx, lnkidx] = Inf
            cappmat[lnkidx, rndidx] = Inf
            if sum(cappconnmat[lnkidx, :]) >= 3
                cappmat[lnkidx, :] .= Inf
                cappmat[:, lnkidx] .= Inf
            end
        end
    end

    # Triangular extraction of the new capp-capp edges.
    connMap = Dict{Tuple{Int,Int}, Int}()
    @inbounds for i in 1:ncapp, j in (i + 1):ncapp
        if cappconnmat[i, j]
            connIdx += 1
            ns = connidxs[i]
            ne = connidxs[j]
            if !(ne in nodes[ns].conn)
                push!(nodes[ns].conn, ne)
                sort!(nodes[ns].conn)
            end
            if !(ns in nodes[ne].conn)
                push!(nodes[ne].conn, ns)
                sort!(nodes[ne].conn)
            end
            push!(conn, gen_conn(ns, ne, NaN, zeros(Int, 0, 3), :capp))
            connMap[(ns, ne)] = connIdx
            connMap[(ne, ns)] = connIdx
        end
    end
    nv = merge(nv, (nconn = connIdx,))

    # Resolve vtcp weights from neighbouring capp weights (BFS).
    vtcp_starts = [c.start for c in conn if c.misc === :vtcp]
    toConnect = Int[]
    for s in vtcp_starts
        for nb in nodes[s].conn
            haskey(connMap, (s, nb)) || continue
            push!(toConnect, connMap[(s, nb)])
        end
    end
    # Add vtcp's own edges into connMap.
    for (i, c) in enumerate(conn)
        if c.misc === :vtcp
            connMap[(c.ends, c.start)] = i
        end
    end

    safety = 100 * length(conn) + 1
    while !isempty(toConnect) && safety > 0
        safety -= 1
        currConn = popfirst!(toConnect)
        cc = conn[currConn]
        isnan(cc.weight) || continue
        startConns = unique([connMap[(cc.start, nb)] for nb in nodes[cc.start].conn
                             if haskey(connMap, (cc.start, nb))])
        endConns   = unique([connMap[(cc.ends, nb)]  for nb in nodes[cc.ends].conn
                             if haskey(connMap, (cc.ends, nb))])
        startConns = filter(x -> x != currConn, startConns)
        endConns   = filter(x -> x != currConn, endConns)
        startW = [conn[k].weight for k in startConns]
        endW   = [conn[k].weight for k in endConns]
        startflag = false
        endflag = false
        w1 = if any(isnan, startW)
            NaN
        elseif length(startW) == 1
            startflag = true
            startW[1]
        elseif !isempty(startW)
            T1 = maximum(startW)^2 - minimum(startW)^2
            T2 = maximum(startW)^2 + minimum(startW)^2
            sqrt(rand(rng) * (T2 - T1) + T1)
        else
            NaN
        end
        w2 = if any(isnan, endW)
            NaN
        elseif length(endW) == 1
            endflag = true
            endW[1]
        elseif !isempty(endW)
            T1 = maximum(endW)^2 - minimum(endW)^2
            T2 = maximum(endW)^2 + minimum(endW)^2
            sqrt(rand(rng) * (T2 - T1) + T1)
        else
            NaN
        end
        connweight = if isnan(w1)
            isnan(w2) ? max(1.0, vesSize[3] + vesSize[4] * randn(rng)) : w2
        else
            isnan(w2) ? w1 : (startflag ? w1 : endflag ? w2 : (w1 + w2) / 2)
        end
        conn[currConn].weight = connweight
        for k in endConns
            isnan(conn[k].weight) && push!(toConnect, k)
        end
        for k in startConns
            isnan(conn[k].weight) && push!(toConnect, k)
        end
    end

    # Fill any leftover NaN weights with a draw from the default capp distribution.
    for c in conn
        if isnan(c.weight)
            c.weight = max(1.0, vesSize[3] + vesSize[4] * randn(rng))
        end
    end

    return nodes, conn, nv
end

# ---------------------------------------------------------------------------
# Render edges into a volume (cscvn → linear-interp simplification)
# ---------------------------------------------------------------------------

"""
    conn_to_vol!(nodes, conn, nv::NamedTuple; idxs=1:length(conn), neur_ves=falses(nv.size))
        -> (neur_ves, conn)

Render the edges indexed by `idxs` into the boolean volume `neur_ves`,
inflating each edge into a tube of radius `conn[i].weight`. Mutates
`conn` to populate `locs`. Ports `connToVol.m`; **the upstream cubic-spline
smoothing (`cscvn`) is replaced by linear interpolation between endpoints**
— the dilation that follows dominates the resulting mask.
"""
function conn_to_vol!(nodes::Vector{VesselNode},
                      conn::Vector{VesselEdge},
                      nv::NamedTuple;
                      idxs=1:length(conn),
                      neur_ves::AbstractArray{Bool,3}=falses(nv.size))
    sz = nv.size
    for i in idxs
        c = conn[i]
        ps = nodes[c.start].pos
        pe = nodes[c.ends].pos
        d = sqrt(sum((ps .- pe).^2))
        numsamp = max(2, Int(ceil(2 * d)))
        ves_loc_set = Set{Tuple{Int,Int,Int}}()
        for k in 1:numsamp
            t = (k - 1) / (numsamp - 1)
            x = Int(ceil((1 - t) * ps[1] + t * pe[1]))
            y = Int(ceil((1 - t) * ps[2] + t * pe[2]))
            z = Int(ceil((1 - t) * ps[3] + t * pe[3]))
            x = clamp(x, 1, sz[1])
            y = clamp(y, 1, sz[2])
            z = clamp(z, 1, sz[3])
            push!(ves_loc_set, (x, y, z))
        end
        ves_loc = collect(ves_loc_set)
        if !isempty(ves_loc)
            mat = zeros(Int, length(ves_loc), 3)
            for (k, t) in enumerate(ves_loc)
                mat[k, 1], mat[k, 2], mat[k, 3] = t
            end
            conn[i].locs = mat
            r = isnan(c.weight) ? 1.0 : max(0.0, c.weight)
            for k in 1:length(ves_loc)
                paint_ball3d!(neur_ves, @view(mat[k, :]), r, true)
            end
        end
    end
    return neur_ves, conn
end

# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

"""
    simulate_blood_vessels(vol_params, vasc_params; rng=Random.default_rng())
        -> (neur_ves, vasc_params)

Generate the in-volume blood-vessel mask. `neur_ves` is a 3-D `BitArray`
whose `true` voxels mark vasculature; the volume is sized
`(vol_sz[1] + 0, vol_sz[2] + 0, vol_sz[3] + vol_depth)` at `vres` voxels/µm
(matching upstream's default `vol_sz + [0 0 vol_depth]` extension).

Ports `simulatebloodvessels.m`. The cubic-spline `cscvn` step is replaced
by linear interpolation; otherwise the algorithm mirrors upstream.
"""
function simulate_blood_vessels(vol_params::VolumeParams,
                                vasc_params::VasculatureParams;
                                rng::AbstractRNG=Random.default_rng())
    vres = vol_params.vres

    # Scaled vasculature parameters.
    vp = (depth_surf     = Int(round(vasc_params.depth_surf * vres)),
          mindists       = vasc_params.vesFreq .* vres ./ 2,
          maxcappdist    = 2 * vasc_params.vesFreq[3] * vres,
          vesSize        = vasc_params.vesSize .* vres,
          vesFreq        = vasc_params.vesFreq,
          sepweight      = vasc_params.sepweight,
          distWeightScale = vasc_params.distWeightScale,
          randWeightScale = vasc_params.randWeightScale,
          distsc         = vasc_params.distsc,
          vesNumScale    = vasc_params.vesNumScale,
          sourceFreq     = vasc_params.sourceFreq)

    # Node-placement parameters, scaled to voxels.
    np_src = vasc_params.node_params
    np = VasculatureNodeParams(
        maxit   = np_src.maxit,
        lensc   = np_src.lensc * vres,
        varsc   = np_src.varsc * vres,
        mindist = np_src.mindist * vres,
        varpos  = np_src.varpos * vres,
        dirvar  = np_src.dirvar,
        branchp = np_src.branchp,
        vesrad  = Int(ceil(np_src.vesrad * vres)),
    )

    # Volume size & node counts.
    vol_sz = Float64[vol_params.vol_sz[1], vol_params.vol_sz[2],
                     vol_params.vol_sz[3] + vol_params.vol_depth]
    sz = (Int(round(vol_sz[1] * vres)),
          Int(round(vol_sz[2] * vres)),
          Int(round(vol_sz[3] * vres)))
    nsource = max(0, Int(round((2 * (vol_sz[1] + vol_sz[2]) / vp.sourceFreq) *
                               abs(1 + vp.vesNumScale * randn(rng)))))
    nvert = max(0, Int(round((vol_sz[1] * vol_sz[2] / vp.vesFreq[2]^2) *
                             abs(1 + vp.vesNumScale * randn(rng)))))
    nsurf = max(0, Int(round((vol_sz[1] * vol_sz[2] / vp.vesFreq[1]^2) *
                             abs(1 + vp.vesNumScale * randn(rng)))))
    ncapp = max(0, Int(round((prod(vol_sz) / vp.vesFreq[3]^3) *
                             abs(1 + vp.vesNumScale * randn(rng)))))

    nv = (size    = sz,
          szum    = (Int(vol_sz[1]), Int(vol_sz[2]), Int(vol_sz[3])),
          nsource = nsource,
          nvert   = nvert,
          nsurf   = nsurf,
          ncapp   = ncapp,
          nnodes  = 0,
          nconn   = 0)

    nodes, nv = grow_major_vessels!(np, vp, nv; rng)

    conn = nodes_to_conn(nodes)
    nv = merge(nv, (nconn = length(conn),))

    # Shift surface vessel z by half of the edge weight (so thick edge nodes don't poke out).
    for c in conn
        for which in (c.start, c.ends)
            t = nodes[which].type
            if t === :edge || t === :surf || t === :sfvt
                nodes[which].pos[3] = min(nodes[which].pos[3] +
                    ceil(c.weight / max(1, length(nodes[which].conn))), sz[3])
            end
        end
    end

    neur_ves = falses(sz)
    conn_to_vol!(nodes, conn, nv; neur_ves)

    nodes, conn, nv = grow_capillaries!(nodes, conn, neur_ves, nv, vp, vres; rng)

    cappidxs = [i for i in 1:length(conn) if isempty(conn[i].locs)]
    conn_to_vol!(nodes, conn, nv; idxs=cappidxs, neur_ves)

    return neur_ves, vasc_params
end
