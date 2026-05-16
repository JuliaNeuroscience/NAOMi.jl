# Spike-train generation.
#
# Derived from `markpointproc.m`, `gen_burst_spike_times.m`,
# `binSpikeTrains.m` in upstream NAOMi-Sim
# (Copyright 2021 Alex Song, Adam Charles, MIT).

export sample_firing_rates, generate_burst_spike_times,
       sample_marked_point_process, bin_spike_trains

"""
    sample_firing_rates([rng,] so::SpikeOptions) -> Vector{Float64}

Draw a length-`so.K` vector of per-neuron firing rates from the distribution
selected by `so.rate_dist`:

- `:gamma` (default): rate ~ Gamma(α = `so.alpha`, scale = `so.rate`), then
  clamped to `[so.rate/10, 10·so.rate]`. With α = 1 this is the exponential
  distribution upstream documents as the default.
- `:uniform` (or anything else): all rates equal to `so.rate`.

Faithful to the preprocessing step in `gen_burst_spike_times.m`.
"""
function sample_firing_rates(rng::AbstractRNG, so::SpikeOptions)
    K = so.K
    base = so.rate
    if so.rate_dist === :gamma
        d = Gamma(so.alpha, base)
        r = rand(rng, d, K)
        clamp!(r, base / 10, base * 10)
        return r
    else  # :uniform and fallback
        return fill(float(base), K)
    end
end

sample_firing_rates(so::SpikeOptions) = sample_firing_rates(Random.default_rng(), so)

"""
    generate_burst_spike_times([rng,] so::SpikeOptions[, rates]) -> Matrix{Int}

Generate a `(K, nt)` binary spike-count matrix following upstream
`gen_burst_spike_times.m`. Each neuron fires as a Poisson process with
per-neuron rate from [`sample_firing_rates`](@ref) (or the supplied `rates`
vector). When `so.burst_mean > 0`, each Poisson event is followed by
`Poisson(so.burst_mean)` additional spikes jittered by `5 + 2·rand()`
sampling intervals (5-sample refractory + uniform 0–2 jitter, as upstream).

Cells of the returned matrix are `0` / `1`; upstream writes `S(k,bin) = 1`
rather than incrementing, so coincident bursts within a single bin do not
accumulate (mirrored here).
"""
function generate_burst_spike_times(rng::AbstractRNG, so::SpikeOptions,
                                    rates::AbstractVector)
    length(rates) == so.K ||
        throw(ArgumentError("rates length $(length(rates)) ≠ so.K = $(so.K)"))
    K, nt = so.K, so.nt
    burst_mean = so.burst_mean
    ref_time = 5.0                       # refractory delay in samples (upstream)
    burst_dist = burst_mean > 0 ? Poisson(burst_mean) : nothing
    S = zeros(Int, K, nt)
    for k in 1:K
        inv_rate = 1.0 / rates[k]        # average inter-spike interval
        T_tot = 0.0
        while T_tot < nt
            t_arr = -inv_rate * log(rand(rng))
            place = min(ceil(Int, T_tot + t_arr), nt + 1)
            if place ≤ nt
                S[k, place] = 1
            end
            T_tot += t_arr
            if burst_dist !== nothing
                num_in_burst = 1 + rand(rng, burst_dist)
                for _ in 2:num_in_burst
                    t_b = ref_time + 2 * rand(rng)
                    place_b = min(ceil(Int, T_tot + t_b), nt + 1)
                    if place_b ≤ nt
                        S[k, place_b] = 1
                    end
                    T_tot += t_b
                end
            end
        end
    end
    return S
end

function generate_burst_spike_times(rng::AbstractRNG, so::SpikeOptions)
    return generate_burst_spike_times(rng, so, sample_firing_rates(rng, so))
end

generate_burst_spike_times(so::SpikeOptions) =
    generate_burst_spike_times(Random.default_rng(), so)

"""
    sample_marked_point_process([rng,]; cif, cifmax=nothing, mkf=nothing,
                                 timemax=Inf, nummax=typemax(Int),
                                 markdim::Integer=0, histlen=Inf)
        -> (times::Vector{Float64}, marks::Matrix{Float64})

Sample a marked point process via Ogata's thinning algorithm (the upstream
`markpointproc.m`).

`cif(t, past_t, past_m)::Real` returns the conditional intensity at time
`t` given past event times `past_t` and marks `past_m`.

`cifmax(t, past_t, past_m) -> (rate_max, interval)` returns an upper bound
on `cif` over the interval `(t, t+interval)`. When `cifmax === nothing` and
the CIF is monotonically nonincreasing between events, the implementation
falls back to evaluating `cif` at `nextfloat(t)` and using
`interval = histlen`, matching the upstream fallback.

`mkf(t, past_t, past_m) -> AbstractVector` produces a `markdim`-length mark
for an event at time `t`; required iff `markdim > 0`.

`histlen` limits how far back the CIF / MKF "see" — older events are
dropped from the `past_t` / `past_m` slices.

Stops when either `timemax` or `nummax` is reached; at least one must be
finite.
"""
function sample_marked_point_process(rng::AbstractRNG;
        cif, cifmax = nothing, mkf = nothing,
        timemax = Inf, nummax::Integer = typemax(Int),
        markdim::Integer = 0, histlen = Inf)
    (isinf(timemax) && nummax == typemax(Int)) &&
        throw(ArgumentError("at least one of timemax or nummax must be finite"))
    markdim ≥ 0 || throw(ArgumentError("markdim must be ≥ 0"))
    markdim == 0 || mkf !== nothing ||
        throw(ArgumentError("mkf must be supplied when markdim > 0"))

    cap = nummax == typemax(Int) ? 100 : Int(nummax)
    evt = zeros(Float64, cap)
    evm = zeros(Float64, cap, markdim)

    _cifmax = cifmax === nothing ?
        (t, past_t, past_m) -> (cif(nextfloat(float(t)), past_t, past_m), histlen) :
        cifmax

    histstart = 1
    t = 0.0
    i = 0
    while i < nummax
        past_t = @view evt[histstart:i]
        past_m = @view evm[histstart:i, :]
        cifmaxt, cifmaxintvl = _cifmax(t, past_t, past_m)
        tstep = -log(rand(rng)) / cifmaxt
        t1 = t + tstep
        if t1 > timemax
            break
        elseif tstep ≥ cifmaxintvl
            t = t + cifmaxintvl
        else
            t = t1
            cift = cif(t, past_t, past_m)
            ratefrac = cift / cifmaxt
            ratefrac > 1 && error("CIF evaluated above CIFMAX (ratefrac=$ratefrac)")
            if rand(rng) < ratefrac
                i += 1
                if i > length(evt)
                    new_cap = 2 * length(evt)
                    evt = vcat(evt, zeros(Float64, new_cap - length(evt)))
                    evm = vcat(evm, zeros(Float64, new_cap - size(evm, 1), markdim))
                end
                evt[i] = t
                if markdim > 0
                    m = mkf(t, past_t, past_m)
                    length(m) == markdim ||
                        error("mkf returned $(length(m)) elements, expected markdim=$markdim")
                    @views evm[i, :] .= m
                end
                while t - histlen > evt[histstart]
                    histstart += 1
                end
            end
        end
    end

    return evt[1:i], evm[1:i, :]
end

sample_marked_point_process(; kwargs...) =
    sample_marked_point_process(Random.default_rng(); kwargs...)

"""
    bin_spike_trains(evt, evm, N_node, dt, T) -> Matrix{Int}

Bin a flat list of marked events into an `N_node × T` count matrix.
`evt[k]` is the (continuous) time of event `k`; `evm[k]` is its integer
mark (the neuron / source index in `1:N_node`); `dt` is the bin width
(same units as `evt`); `T` is the number of output bins.

Bin index for event `k` is `ceil(evt[k]/dt)`; counts accumulate, so
coincident events within a bin sum.
"""
function bin_spike_trains(evt::AbstractVector{<:Real},
                          evm::AbstractVector{<:Integer},
                          N_node::Integer, dt::Real, T::Integer)
    length(evt) == length(evm) ||
        throw(ArgumentError("evt and evm must have the same length"))
    if !isempty(evm) && maximum(evm) > N_node
        throw(ArgumentError("largest mark exceeds N_node"))
    end
    if !isempty(evt) && ceil(Int, maximum(evt) / dt) > T
        throw(ArgumentError("latest event ($(maximum(evt))) exceeds T·dt = $(T*dt)"))
    end
    S = zeros(Int, N_node, T)
    for k in eachindex(evt)
        bin = ceil(Int, evt[k] / dt)
        S[evm[k], bin] += 1
    end
    return S
end
