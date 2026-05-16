# Calcium dynamics and fluorescence transduction.
#
# Derived from `make_calcium_impulse.m`, `calcium_dynamics.m`, and
# `mk_doub_exp_ker.m` in upstream NAOMi-Sim
# (Copyright 2021 Alex Song, Adam Charles, MIT). The streaming variants
# (`genNextCalciumDynamics.m`, `genNextSpikeTimepoint.m`,
# `generateNextTimePoint.m`) are *not* ported here — they exist solely to
# support the LowRAM script variant, which is explicitly out of scope (see
# ANALYSIS_PLAN.md "Out of scope").

export make_doub_exp_kernel, make_calcium_impulse, calcium_dynamics,
       fluorescence

"""
    make_doub_exp_kernel(t_on, t_off, A, dt) -> Vector{Float64}

Causal double-exponential kernel sampled at intervals `dt`, evaluating

    h(t) = A · (1 - exp(-t_on · t)) · exp(-t_off · t)

from `t = 0` up to where the kernel decays below `1e-3` of its peak. Note
that `t_on` and `t_off` are *rates* (units of `1/time`), matching the
upstream MATLAB convention despite the misleading names; with the
`CalciumParams` defaults (`t_on = 0.8535`, `t_off = 98.6173` for GCaMP6f)
the peak occurs near 50 ms.

Implements the default ('mult' / `otherwise`) branch of upstream
`mk_doub_exp_ker.m`; the `'plus'` and `'min'` branches are deferred
(they are not exercised by the standard pipeline).
"""
function make_doub_exp_kernel(t_on::Real, t_off::Real, A::Real, dt::Real)
    t_on > 0 && t_off > 0 ||
        throw(ArgumentError("t_on and t_off must be positive"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
    loc_max = log((t_off + t_on) / t_off) / t_on
    max_val = A * (1 - exp(-t_on * loc_max)) * exp(-t_off * loc_max)
    t_max = -log(max_val * 1e-3 / A) / t_off
    ts = 0:dt:(t_max + dt)
    return [A * (1 - exp(-t_on * t)) * exp(-t_off * t) for t in ts]
end

"""
    make_calcium_impulse(ca_scale; dt = 1/30) -> Vector{Float64}

Impulse response of the all-pole AR system whose poles are
`exp.(-ca_scale)`, sampled for `ceil(10/dt)` steps. Ports upstream
`make_calcium_impulse.m`, which uses MATLAB's `arima` / `impulse`; the
equivalent here is a direct difference-equation evaluation, avoiding any
control-theory dependency. `ca_scale` may be any length; the standard
two-element case `[t_on, t_off]` yields a two-pole low-pass response.
"""
function make_calcium_impulse(ca_scale::AbstractVector{<:Real};
                              dt::Real = 1 / 30)
    maximum(ca_scale) > 0 ||
        throw(ArgumentError("ca_scale must contain at least one positive value"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
    poles = exp.(-float.(ca_scale))
    coeffs = _poly_from_roots(poles)        # [1, a1, a2, …, aN]
    N = length(coeffs) - 1
    L = ceil(Int, 10 / dt)
    h = zeros(Float64, L)
    h[1] = 1.0                              # impulse at sample 1
    for n in 2:L
        s = 0.0
        for k in 1:N
            n - k ≥ 1 && (s -= coeffs[k + 1] * h[n - k])
        end
        h[n] = s
    end
    return h
end

# Expand prod_i (z - r_i) into descending-power coefficients [1, a1, ..., aN]
# so that the AR recursion is y[n] = -a1·y[n-1] - … - aN·y[n-N] + x[n].
function _poly_from_roots(roots::AbstractVector)
    a = [1.0]
    for r in roots
        b = zeros(Float64, length(a) + 1)
        for i in eachindex(a)
            b[i]     += a[i]
            b[i + 1] -= r * a[i]
        end
        a = b
    end
    return a
end

"""
    calcium_dynamics(S, cp::CalciumParams; over_samp = 1, ext_mult = 1)
        -> (; CB, C, F)

Simulate the protein ↔ Ca²⁺ binding model driven by the spike-count matrix
`S` (`K × nt`) and return total calcium `C`, bound-indicator `CB`, and
fluorescence `F` (each `K × nt`).

The dynamics flavour is selected by `cp.sat_type`:

- `:double` (upstream default) — two-binding-site model. The bound state is
  the sum of fast and slow components.
- `:single` — single-binding-site model. Before passing `CB` to
  [`fluorescence`](@ref) it is rescaled to `CB + ca_rest + (b/a)·CB/(ind_con − CB)`
  to recover the bound-Ca interpretation the Hill equation expects.
- `:Ca_DE` — single-pool ODE (no explicit binding sites) followed by
  convolution with a double-exponential kernel built by
  [`make_doub_exp_kernel`](@ref). Because that branch already smooths via
  convolution, no extra resampling is applied.

`over_samp > 1` repeats each spike `over_samp` times before integrating
(upstream's oversampling shortcut); the returned arrays are decimated back
to the original column count. `ext_mult` scales the extrusion rate.

When `cp.ca_sat ∈ [0, 1)`, total calcium `C` is capped at
`ca_dis · ca_sat / (1 - ca_sat)` per upstream.
"""
function calcium_dynamics(S::AbstractMatrix{<:Real}, cp::CalciumParams;
                          over_samp::Integer = 1, ext_mult::Real = 1)
    over_samp ≥ 1 || throw(ArgumentError("over_samp must be ≥ 1"))
    K, nt_in = size(S)

    if over_samp > 1
        Sx = zeros(Float64, K, nt_in * over_samp)
        for j in 1:nt_in
            Sx[:, (j - 1) * over_samp + 1] = S[:, j]
        end
    else
        Sx = Float64.(S)
    end
    nt = size(Sx, 2)

    ext_rate = ext_mult * cp.ext_rate
    ca_bind  = cp.ca_bind
    ca_rest  = cp.ca_rest
    ind_con  = cp.ind_con
    ca_dis   = cp.ca_dis
    ca_sat   = cp.ca_sat
    dt       = cp.dt

    sat_cap = (0 ≤ ca_sat < 1) ? ca_dis * ca_sat / (1 - ca_sat) : Inf

    C = zeros(Float64, K, nt)
    @views C[:, 1] .= max.(ca_rest, Sx[:, 1])

    local CB_out::Matrix{Float64}
    needs_decimation = true

    if cp.sat_type === :single
        a, b = cp.a_bind, cp.a_ubind
        CB = zeros(Float64, K, nt)
        @inbounds for j in 2:nt
            for k in 1:K
                Cprev = C[k, j - 1]
                CBprev = CB[k, j - 1]
                denom = 1 + ca_bind + (ind_con * ca_dis) / (Cprev + ca_dis)^2
                Cnew = Cprev + dt * b * CBprev +
                       (-dt * ext_rate * (Cprev - CBprev - ca_rest) + Sx[k, j]) / denom
                C[k, j] = min(Cnew, sat_cap)
                CB[k, j] = CBprev + dt * (-b * CBprev +
                           a * (Cprev - CBprev) * (ind_con - CBprev))
            end
        end
        CB_out = CB
    elseif cp.sat_type === :Ca_DE
        @inbounds for j in 2:nt
            for k in 1:K
                Cprev = C[k, j - 1]
                denom = 1 + ca_bind + (ind_con * ca_dis) / (Cprev + ca_dis)^2
                Cnew = Cprev + (-dt * ext_rate * (Cprev - ca_rest) + Sx[k, j]) / denom
                C[k, j] = min(Cnew, sat_cap)
            end
        end
        h = make_doub_exp_kernel(cp.t_on, cp.t_off, cp.ca_amp, dt)
        CB = similar(C)
        for k in 1:K
            CB[k, :] .= _conv_full_decimate(view(C, k, :) .- ca_rest, h,
                                            over_samp, nt_in) .+ ca_rest
        end
        # Decimate C too; Ca_DE branch does its own resampling here.
        if over_samp > 1
            C = C[:, 1:over_samp:end]
        end
        CB_out = CB
        needs_decimation = false
    elseif cp.sat_type === :double
        a = _aspair(cp.a_bind)
        b = _aspair(cp.a_ubind)
        CB1 = zeros(Float64, K, nt)
        CB2 = zeros(Float64, K, nt)
        @inbounds for j in 2:nt
            for k in 1:K
                Cprev   = C[k, j - 1]
                CB1prev = CB1[k, j - 1]
                CB2prev = CB2[k, j - 1]
                denom = 1 + ca_bind + (ind_con * ca_dis) / (Cprev + ca_dis)^2
                Cnew = Cprev + dt * (b[1] * CB1prev + b[2] * CB2prev) +
                       (-dt * ext_rate * (Cprev - CB1prev - CB2prev - ca_rest) +
                        Sx[k, j]) / denom
                C[k, j] = min(Cnew, sat_cap)
                free = Cprev - CB1prev - CB2prev
                bound_room = ind_con - CB1prev - CB2prev
                CB1[k, j] = CB1prev + dt * (-b[1] * CB1prev + a[1] * free * bound_room)
                CB2[k, j] = CB2prev + dt * (-b[2] * CB2prev + a[2] * free * bound_room)
            end
        end
        CB_out = CB1 .+ CB2
    else
        throw(ArgumentError("unknown sat_type: $(cp.sat_type)"))
    end

    if needs_decimation && over_samp > 1
        C = C[:, 1:over_samp:end]
        CB_out = CB_out[:, 1:over_samp:end]
    end

    CB_for_F = cp.sat_type === :single ?
        CB_out .+ ca_rest .+ (cp.a_ubind / cp.a_bind) .* CB_out ./ (ind_con .- CB_out) :
        CB_out
    F = fluorescence(CB_for_F, cp.prot)
    return (; CB = CB_out, C = C, F = F)
end

_aspair(x::Real) = (float(x), float(x))
function _aspair(x::AbstractVector{<:Real})
    n = length(x)
    n == 1 && return (float(x[1]), float(x[1]))
    return (float(x[1]), float(x[2]))
end

# Causal "full" convolution of `signal` with `kernel`, decimated by
# `over_samp` and truncated to `nt_in` columns. Hand-rolled to avoid a
# DSP.jl dependency.
function _conv_full_decimate(signal::AbstractVector{<:Real},
                             kernel::AbstractVector{<:Real},
                             over_samp::Integer, nt_in::Integer)
    Ls, Lk = length(signal), length(kernel)
    Lf = Ls + Lk - 1
    out_len = min(nt_in, cld(Lf, over_samp))
    out = zeros(Float64, out_len)
    for j in 1:out_len
        n = (j - 1) * over_samp + 1   # 1-based index in the full conv
        s = 0.0
        kmin = max(1, n - Ls + 1)
        kmax = min(Lk, n)
        for k in kmin:kmax
            s += kernel[k] * signal[n - k + 1]
        end
        out[j] = s
    end
    return out
end

"""
    fluorescence(CB, prot::Symbol) -> Array

Hill-equation fluorescence transduction. `CB` is bound-indicator
concentration (mol/L); the return has the same shape as `CB` and
represents the calibrated fluorescence ``F = F_0 (1 + dF/F)`` per
Badura et al. 2014 (GCaMP6 / GCaMP3 / OGB-1, RS06, RS09) and
Dana et al. 2019 (jGCaMP7 family, GCaMP6s).

Accepted protein symbols (case-sensitive aliases mirror upstream
`switch lower(prot_type)`):

- `:GCaMP6`, `:GCaMP6f`
- `:GCaMP6s`
- `:GCaMP3`
- `:OGB1`, `:OGB_1`
- `:GCaMP6_RS06`, `:GCaMP6rs06`
- `:GCaMP6_RS09`, `:GCaMP6rs09`
- `:jGCaMP7f`, `:jGCaMP7s`, `:jGCaMP7b`, `:jGCaMP7c`

Unknown symbols emit a one-line warning and fall back to the `:GCaMP6f`
parameters, mirroring the upstream `otherwise` branch.
"""
function fluorescence(CB::AbstractArray{<:Real}, prot::Symbol)
    F0, Fmax, K_d, h = _fluor_params(prot)
    return @. F0 + F0 * Fmax / (1 + (K_d / CB)^h)
end

function _fluor_params(prot::Symbol)
    if prot === :GCaMP6 || prot === :GCaMP6f
        return (1.0, 25.2, 290e-9, 2.7)
    elseif prot === :GCaMP6s
        return (1.0, 27.2, 147e-9, 2.45)
    elseif prot === :GCaMP3
        return (2.0, 12.0, 287e-9, 2.52)
    elseif prot === :OGB1 || prot === :OGB_1
        return (1.0, 14.0, 250e-9, 1.0)
    elseif prot === :GCaMP6_RS09 || prot === :GCaMP6rs09
        return (1.4, 25.0, 520e-9, 3.2)
    elseif prot === :GCaMP6_RS06 || prot === :GCaMP6rs06
        return (1.2, 15.0, 320e-9, 3.0)
    elseif prot === :jGCaMP7f
        return (1.0, 30.2, 174e-9, 2.3)
    elseif prot === :jGCaMP7s
        return (1.0, 40.4, 68e-9, 2.49)
    elseif prot === :jGCaMP7b
        return (1.0, 22.1, 82e-9, 3.06)
    elseif prot === :jGCaMP7c
        return (1.0, 145.6, 298e-9, 2.44)
    else
        @warn "Unknown protein $(prot); defaulting to :GCaMP6f"
        return (1.0, 25.2, 290e-9, 2.7)
    end
end
