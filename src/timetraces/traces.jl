# Top-level time-trace orchestration.
#
# Derived from `generateTimeTraces.m`, `genCorrelatedSpikeTrains2.m`,
# `expression_variation.m`, and `sampSmallWorldMat.m` in upstream NAOMi-Sim
# (Copyright 2021 Alex Song, Adam Charles, MIT). The continuous-time
# (`discrete_flag = false`) path through `markpointproc` is *not* ported;
# the standard pipeline always uses the discrete approximation.

using StatsBase: sample

export samp_small_world_mat, expression_variation,
       gen_correlated_spike_trains, generate_time_traces

"""
    samp_small_world_mat([rng,] N_node, K_conn, beta;
                          rand_opt = 0, self_ex = 4, n_locs = nothing)
        -> Matrix{Float64}

Watts–Strogatz-style small-world adjacency matrix. `N_node` may be an
integer (no background nodes) or a 2-tuple `(K, N_bg)` (background nodes
appended; they receive a fully connected fan-in from the soma block).

`K_conn` is the per-node lattice connection count (must be even). `beta`
is the rewiring probability. `rand_opt` ∈ [0, 1] mixes in
uniform-random weights on `[0.1, 1.0]`; `self_ex` is added to the
diagonal to drive bursting.

`n_locs` (if supplied) is a `K × d` location matrix; when present the
initial `K_conn` connections of each node are its `K_conn` nearest
spatial neighbours rather than a Toeplitz lattice.
"""
function samp_small_world_mat(rng::AbstractRNG, N_node, K_conn::Integer, beta;
                              rand_opt = 0.0, self_ex = 4.0, n_locs = nothing)
    K, N_bg = N_node isa Tuple || N_node isa AbstractVector ?
        (Int(N_node[1]), Int(N_node[2])) : (Int(N_node), 0)
    iseven(K_conn) || throw(ArgumentError("K_conn must be even"))
    K ≥ K_conn || throw(ArgumentError("K_conn ($(K_conn)) exceeds N_node ($K)"))

    use_locs = n_locs !== nothing && size(n_locs, 1) == K

    adj = if use_locs
        # Pairwise Euclidean distance (K×K)
        D = zeros(Float64, K, K)
        for d in 1:size(n_locs, 2), j in 1:K, i in 1:K
            D[i, j] += (n_locs[i, d] - n_locs[j, d])^2
        end
        D .= sqrt.(D)
        out = zeros(Float64, K, K)
        for i in 1:K
            order = sortperm(@view D[i, :])
            for j in order[1:K_conn]
                out[i, j] = 1.0
            end
        end
        out
    else
        # 1-D lattice (matches upstream's `toeplitz` init — not a ring)
        first_row = [j ≤ K_conn ÷ 2 ? 1.0 : 0.0 for j in 1:K]
        [first_row[abs(i - j) + 1] for i in 1:K, j in 1:K]
    end

    for i in 1:K
        ones_pos  = findall(==(1.0), adj[i, :])
        zeros_pos = findall(==(0.0), adj[i, :])
        n_switch  = count(rand(rng, length(ones_pos)) .< beta)
        n_switch == 0 && continue
        n_switch  = min(n_switch, length(zeros_pos), length(ones_pos))
        new_idx   = sample(rng, 1:length(zeros_pos), n_switch; replace = false)
        # Mirror upstream: switch off the first `n_switch` "on" connections
        # in row order (upstream draws the per-slot mask uniformly, so the
        # count alone determines the rewire effect).
        for k in 1:n_switch
            adj[i, ones_pos[k]] = 0.0
            adj[i, zeros_pos[new_idx[k]]] = 1.0
        end
    end

    if N_bg > 0
        N_total = K + N_bg
        new_adj = zeros(Float64, N_total, N_total)
        new_adj[1:K, 1:K] .= adj
        new_adj[K + 1:end, 1:K] .= 1.0           # bg ← all somas
        for i in 1:N_bg
            new_adj[K + i, K + i] = 1.0
        end
        adj = new_adj
    end

    for i in axes(adj, 1)
        adj[i, i] += self_ex
    end
    if rand_opt != 0
        adj .*= (1 - rand_opt) .+ rand_opt .* (0.1 .+ 0.9 .* rand(rng, size(adj)...))
    end
    return adj
end

samp_small_world_mat(N_node, K_conn::Integer, beta; kwargs...) =
    samp_small_world_mat(Random.default_rng(), N_node, K_conn, beta; kwargs...)

"""
    expression_variation([rng,] N::Integer, p_off, min_mod) -> Vector{Float64}

Per-cell multiplicative expression factor. With probability `p_off`, the
cell is silenced (factor `0`); otherwise the factor is drawn from
either a uniform on `[min_mod, 1]` (scalar `min_mod`) or a
`Gamma(min_mod[2], min_mod[1])` (2-element `min_mod`, matching
upstream's `gamrnd(min_mod(2), min_mod(1))`).
"""
function expression_variation(rng::AbstractRNG, N::Integer, p_off::Real, min_mod)
    0 ≤ p_off ≤ 1 || throw(ArgumentError("p_off must lie in [0, 1]"))
    if min_mod isa Real
        (0 ≤ min_mod ≤ 1) || throw(ArgumentError("min_mod must lie in [0, 1]"))
        x = min_mod .+ (1 - min_mod) .* rand(rng, N)
    elseif length(min_mod) == 1
        m = min_mod[1]
        (0 ≤ m ≤ 1) || throw(ArgumentError("min_mod must lie in [0, 1]"))
        x = m .+ (1 - m) .* rand(rng, N)
    else
        x = rand(rng, Gamma(min_mod[2], min_mod[1]), N)
    end
    return x .* (p_off .< rand(rng, N))
end

expression_variation(N::Integer, p_off::Real, min_mod) =
    expression_variation(Random.default_rng(), N, p_off, min_mod)

"""
    gen_correlated_spike_trains([rng,] spike_opts::SpikeOptions;
                                 n_locs = nothing)
        -> (; soma, bg, A, MU, B)

Run the upstream-style discrete Hawkes process for `spike_opts.nt`
time-steps at `spike_opts.dt`, returning soma and background spike
matrices (`K × nt` and `N_bg × nt`, integer counts). Connectivity is
seeded by [`samp_small_world_mat`](@ref). The continuous-time
(`discrete_flag = false`) branch is deferred — the standard pipeline
uses the discrete approximation by default.

The excitation matrix `A`, baseline rates `MU`, and self-inhibition
constants `B` are returned alongside for inspection.
"""
function gen_correlated_spike_trains(rng::AbstractRNG, so::SpikeOptions;
                                     n_locs = nothing)
    K, N_bg, dt, nt = so.K, so.N_bg, so.dt, so.nt
    N_tot = K + N_bg
    ascale, bscale = 4.0, 2.0
    A = samp_small_world_mat(rng, (K, N_bg), 10, 0.3;
                             rand_opt = 0.9, self_ex = so.burst_mean,
                             n_locs = n_locs)
    A .= ascale .* A ./ (sum(A) / size(A, 1))   # mean of row-sums
    MU = rand(rng, Gamma(1.0, so.rate), N_tot)
    B  = rand(rng, Gamma(3.0, bscale), N_tot)
    for i in 1:N_tot
        A[i, i] = so.selfact * B[i]
    end

    # Upstream sizes extSc/inbSc to K (a latent bug when N_bg > 0); we
    # size to N_tot so the function is correct for every N_bg.
    extSc = max.(0.3, 1.0 .+ 0.3 .* randn(rng, N_tot))
    inbSc = extSc ./ 2
    alpha = 3.0

    yt = fill(5.0, N_tot)
    zt = zeros(Float64, N_tot)
    S = zeros(Int, N_tot, nt)
    for tt in 1:nt
        z = zt .- yt .+ 1.0
        rect = log.(1.0 .+ exp.(alpha .* z))
        thresh = 1.0 .- exp.(.-rect .* MU .* dt)
        xt = rand(rng, N_tot) .< thresh
        zt = exp.(.-extSc .* dt) .* zt .+ A * xt
        yt = exp.(.-inbSc .* dt) .* yt .+ B .* xt
        for k in 1:N_tot
            xt[k] && (S[k, tt] = 1)
        end
    end
    soma = S[1:K, :]
    bg   = N_bg > 0 ? S[K + 1:end, :] : zeros(Int, 0, nt)
    return (; soma, bg, A, MU, B)
end

gen_correlated_spike_trains(so::SpikeOptions; kwargs...) =
    gen_correlated_spike_trains(Random.default_rng(), so; kwargs...)

"""
    generate_time_traces([rng,] spike_opts::SpikeOptions;
                          cal_params = CalciumParams(spike_opts.prot),
                          n_locs = nothing, mod_vals = nothing,
                          S_times = nothing)
        -> (; soma, dend, bg, spikes, mod_vals)

End-to-end time-trace generator. Returns fluorescence traces for soma,
optionally dendrite (`spike_opts.dendflag`) and axon/background
(`spike_opts.axonflag` or `spike_opts.N_bg > 0`), plus the underlying
spike matrix and the per-cell expression-modulation factors.

Internal simulation runs at 100 Hz (matching upstream); when
`spike_opts.dt != 1/100`, the output is linearly resampled to the
requested rate. Upstream uses a polyphase resampler; linear interp is
a deliberate simplification — adequate for downstream chunks.
"""
function generate_time_traces(rng::AbstractRNG, so::SpikeOptions;
                              cal_params::CalciumParams = CalciumParams(so.prot),
                              n_locs = nothing, mod_vals = nothing,
                              S_times = nothing)
    dt_user = so.dt
    nt_user = so.nt
    so_sim = deepcopy(so)
    so_sim.dt = 1 / 100
    so_sim.nt = ceil(Int, nt_user * 100 * dt_user)
    l_buff = 500

    # 1. Spike generation
    if S_times === nothing
        if so.smod_flag === :hawkes
            res = gen_correlated_spike_trains(rng, so_sim; n_locs = n_locs)
            S_soma = res.soma
            S_bg_in = res.bg
        else
            S_soma = generate_burst_spike_times(rng, so_sim)
            S_bg_in = zeros(Int, 0, so_sim.nt)
        end
    else
        size(S_times, 1) == so.K ||
            throw(ArgumentError("S_times rows ($(size(S_times,1))) ≠ K = $(so.K)"))
        S_soma = Int.(S_times)
        S_bg_in = zeros(Int, 0, size(S_soma, 2))
    end
    spikes_soma = S_soma

    # 2. Buffer / rescale into mol/L
    S_soma_buf = _buffer_and_scale(S_soma, l_buff, so.dyn_type, so.mu, so.sig, rng)

    cp = cal_params
    cp.dt = 1 / 100
    prot = cp.prot

    # 3. Per-compartment traces
    soma_trace = _traces_for_compartment(S_soma_buf, cp, so.dyn_type, prot,
                                         :soma, l_buff, rng)
    dend_trace = so.dendflag ?
        _traces_for_compartment(S_soma_buf, cp, so.dyn_type, prot,
                                :dend, l_buff, rng) : nothing

    bg_trace = nothing
    if so.N_bg > 0
        S_bg = if so.smod_flag === :hawkes
            S_bg_in
        else
            so_bg = deepcopy(so_sim)
            so_bg.K = so.N_bg
            so_bg.rate = 0.25
            so_bg.sig = 0.2
            generate_burst_spike_times(rng, so_bg)
        end
        S_bg_buf = _buffer_and_scale(S_bg, l_buff, so.dyn_type, so.mu, 0.2, rng)
        bg_trace = _traces_for_compartment(S_bg_buf, cp, so.dyn_type, prot,
                                           :bg, l_buff, rng)
    end

    # 4. Resample to user dt
    soma_out = _resample_to_user(soma_trace, dt_user, nt_user)
    dend_out = dend_trace === nothing ? nothing :
               _resample_to_user(dend_trace, dt_user, nt_user)
    bg_out   = bg_trace   === nothing ? nothing :
               _resample_to_user(bg_trace,   dt_user, nt_user)

    # 5. Expression modulation
    K_eff = size(soma_out, 1)
    mod_vals === nothing &&
        (mod_vals = expression_variation(rng, K_eff, so.p_off, so.min_mod))
    length(mod_vals) == K_eff ||
        throw(ArgumentError("mod_vals length ($(length(mod_vals))) ≠ K = $K_eff"))
    soma_out .= soma_out .* mod_vals
    dend_out === nothing || (dend_out .= dend_out .* mod_vals)
    bg_out !== nothing && size(bg_out, 1) == K_eff && (bg_out .= bg_out .* mod_vals)

    return (; soma = soma_out, dend = dend_out, bg = bg_out,
              spikes = spikes_soma, mod_vals)
end

generate_time_traces(so::SpikeOptions; kwargs...) =
    generate_time_traces(Random.default_rng(), so; kwargs...)

# --- internal helpers ----------------------------------------------------

function _buffer_and_scale(S::AbstractMatrix, l_buff::Integer, dyn_type::Symbol,
                           mu::Real, sig::Real, rng::AbstractRNG)
    if dyn_type === :AR1 || dyn_type === :AR2
        out = Float64.(S)
        for i in eachindex(out)
            if out[i] == 1
                out[i] = (1 + rand(rng)) * exp(mu + sig * randn(rng))
            end
        end
        return out                                # AR branches: no l_buff prepend
    else
        K = size(S, 1)
        out = hcat(zeros(Float64, K, l_buff), Float64.(S))
        out .*= 7.6e-6
        return out
    end
end

function _ext_rate_override(dyn_type::Symbol, compartment::Symbol)
    if dyn_type === :single
        compartment === :soma && return 800.0
        compartment === :dend && return 400.0
        compartment === :bg   && return 1600.0
    elseif dyn_type === :double
        compartment === :soma && return 800.0
        compartment === :dend && return 1000.0
        compartment === :bg   && return 200.0
    elseif dyn_type === :Ca_DE
        compartment === :bg && return 2800.0
    end
    return nothing
end

function _ext_mult_override(dyn_type::Symbol, compartment::Symbol)
    if dyn_type === :Ca_DE
        compartment === :dend && return 0.5
        compartment === :bg   && return 0.25
    end
    return 1.0
end

function _traces_for_compartment(S_buf::AbstractMatrix, cp::CalciumParams,
                                 dyn_type::Symbol, prot::Symbol,
                                 compartment::Symbol, l_buff::Integer,
                                 rng::AbstractRNG)
    if dyn_type === :AR1 || dyn_type === :AR2
        # AR1 uses one pole; AR2 uses two. Upstream uses `make_calcium_impulse(0.9, 1)` for soma and (0.8, 1) for dend/bg.
        pole = compartment === :soma ? 0.9 : 0.8
        ca_scale = dyn_type === :AR1 ? [pole] : [pole, pole]
        h = make_calcium_impulse(ca_scale; dt = cp.dt)
        h ./= maximum(h)
        h .*= 0.5
        K = size(S_buf, 1)
        b_cell = abs.(1 .+ 0.1 .* randn(rng, K))
        scale = 0.5 + 0.5 * rand(rng)
        convolved = _conv_same(S_buf, h)
        return 2.5 .* convolved .* (b_cell .* scale) .+ b_cell
    else
        cp_local = deepcopy(cp)
        cp_local.sat_type = dyn_type === :double ? :double :
                            dyn_type === :single ? :single : :Ca_DE
        ext_override = _ext_rate_override(dyn_type, compartment)
        ext_mult     = _ext_mult_override(dyn_type, compartment)
        ext_override !== nothing && (cp_local.ext_rate = ext_override)
        # :single / :double rescale by min(S>0) to keep concentrations in
        # the expected range (upstream behaviour).
        S_in = S_buf
        if dyn_type === :single || dyn_type === :double
            pos = filter(>(0), S_buf)
            if !isempty(pos)
                S_in = (7.6e-6) .* S_buf ./ minimum(pos)
            end
        end
        res = calcium_dynamics(S_in, cp_local; ext_mult)
        return res.F[:, l_buff + 1:end]
    end
end

function _conv_same(S::AbstractMatrix, h::AbstractVector)
    K, T = size(S)
    Lh = length(h)
    out = zeros(Float64, K, T)
    offset = Lh ÷ 2
    for k in 1:K, t in 1:T
        s = 0.0
        for j in 1:Lh
            idx = t + offset - j + 1
            (1 ≤ idx ≤ T) && (s += h[j] * S[k, idx])
        end
        out[k, t] = s
    end
    return out
end

function _resample_to_user(X::AbstractMatrix, dt_user::Real, nt_user::Integer)
    K, T = size(X)
    if isapprox(dt_user, 1 / 100; atol = 1e-12)
        return Float64.(X[:, 1:min(T, nt_user)])
    end
    out = zeros(Float64, K, nt_user)
    for j in 1:nt_user
        x = (j - 1) * dt_user * 100 + 1
        x = clamp(x, 1.0, Float64(T))
        i0 = floor(Int, x); i1 = min(i0 + 1, T)
        frac = x - i0
        for k in 1:K
            out[k, j] = (1 - frac) * X[k, i0] + frac * X[k, i1]
        end
    end
    return out
end
