# Parameter structs for the NAOMi simulation pipeline.
#
# Derived from the upstream `check_*_params.m` / `check_*_opts.m` files of
# NAOMi-Sim (Copyright 2021 Alex Song, Adam Charles, MIT). Field names are
# preserved verbatim from upstream so they map 1:1 to upstream documentation
# and example scripts.

export VolumeParams, NeuronParams, VasculatureParams, VasculatureNodeParams,
       DendriteParams, AxonParams, BackgroundParams,
       SpikeOptions, CalciumParams,
       NoiseParams, PSFParams, PSFFastMask, ScanParams, TPMParams,
       finalize!

# ---------------------------------------------------------------------------
# Volume
# ---------------------------------------------------------------------------

"""
    VolumeParams(; kwargs...)

Parameters for neural-volume generation. `vol_sz` is in microns; `vres` is
voxels-per-micron sampling resolution. `N_neur == 0` (the default) means
"derive from `neur_density`" — call [`finalize!`](@ref) to fill it in.

Field defaults follow `check_vol_params.m` from upstream NAOMi-Sim.
"""
Base.@kwdef mutable struct VolumeParams
    vol_sz::Vector{Int}      = [100, 100, 50]      # Volume size (µm)
    min_dist::Float64        = 16.0                # Minimum inter-neuron distance (µm)
    vres::Float64            = 2.0                 # Voxels per µm
    N_bg::Int                = 1_000_000           # Number of background processes
    vol_depth::Float64       = 200.0               # Depth of volume centre under brain surface (µm)
    dendrite_tau::Float64    = 5.0                 # Dendrite decay strength
    verbose::Int             = 1                   # 0 silent, 1 brief, 2 detailed
    N_neur::Int              = 0                   # 0 → derive from neur_density
    neur_density::Float64    = 1e5                 # Cells per mm³
    N_den::Float64           = 0.0                 # 0 → derive from AD_density
    AD_density::Float64      = 2e3                 # Apical dendrites per mm²
end

# ---------------------------------------------------------------------------
# Neurons
# ---------------------------------------------------------------------------

"""
    NeuronParams(; kwargs...)

Parameters for individual neuron (soma + nucleus) generation. Defaults follow
`check_neur_params.m`.
"""
Base.@kwdef mutable struct NeuronParams
    n_samps::Int             = 200
    l_scale::Float64         = 90.0
    p_scale::Float64         = 1000.0
    avg_rad::Float64         = 5.9
    nuc_rad::Vector{Float64} = [5.65, 2.5]
    max_ang::Float64         = 20.0
    plot_opt::Bool           = false
    dendrite_tau::Float64    = 50.0
    nuc_fluorsc::Float64     = 0.0
    min_thic::Vector{Float64}   = [0.4, 0.4]
    eccen::Vector{Float64}      = [0.35, 0.35, 0.5]
    exts::Vector{Float64}       = [0.75, 1.7]
    nexts::Vector{Float64}      = [0.5, 1.0]
    neur_type::Symbol           = :pyr
    fluor_dist::Vector{Float64} = [1.0, 0.2]
end

# ---------------------------------------------------------------------------
# Vasculature
# ---------------------------------------------------------------------------

"""
    VasculatureNodeParams(; kwargs...)

Vasculature node-placement parameters. Defaults follow the `node_params`
sub-struct in `check_vasc_params.m`.
"""
Base.@kwdef mutable struct VasculatureNodeParams
    maxit::Int       = 25
    lensc::Float64   = 50.0
    varsc::Float64   = 15.0
    mindist::Float64 = 10.0
    varpos::Float64  = 5.0
    dirvar::Float64  = π / 8
    branchp::Float64 = 0.02
    vesrad::Float64  = 25.0
end

"""
    VasculatureParams(; kwargs...)

Parameters for vasculature generation. Defaults follow `check_vasc_params.m`.
"""
Base.@kwdef mutable struct VasculatureParams
    flag::Int                          = 1
    ves_shift::Vector{Float64}         = [5.0, 15.0, 5.0]
    depth_vasc::Float64                = 200.0
    depth_surf::Float64                = 15.0
    distWeightScale::Float64           = 2.0
    randWeightScale::Float64           = 0.1
    cappAmpScale::Float64              = 0.5
    cappAmpZscale::Float64             = 0.5
    vesSize::Vector{Float64}           = [15.0, 9.0, 2.0]
    vesFreq::Vector{Float64}           = [125.0, 200.0, 50.0]
    sourceFreq::Float64                = 1000.0
    vesNumScale::Float64               = 0.2
    sepweight::Float64                 = 0.75
    distsc::Float64                    = 4.0
    node_params::VasculatureNodeParams = VasculatureNodeParams()
end

# ---------------------------------------------------------------------------
# Dendrites
# ---------------------------------------------------------------------------

"""
    DendriteParams(; kwargs...)

Parameters for dendrite generation. Defaults follow `check_dend_params.m`.

`dtParams`, `atParams`, `atParams2` are upstream-style 5-element vectors
documenting dendritic-tree branching geometry; see upstream docstring for
the per-element meaning.
"""
Base.@kwdef mutable struct DendriteParams
    dtParams::Vector{Float64}    = [40.0, 150.0, 50.0, 1.0, 10.0]
    atParams::Vector{Float64}    = [6.0, 5.0, 5.0, 5.0, 1.0]
    atParams2::Vector{Float64}   = [1.0, 5.0, 5.0, 5.0, 4.0]
    dweight::Float64             = 10.0
    bweight::Float64             = 5.0
    thicknessScale::Float64      = 0.5
    weightScale::Vector{Float64} = [150.0, 1.0, 0.8]
    dims::Vector{Int}            = [60, 60, 60]
    dimsSS::Vector{Int}          = [5, 5, 5]
    rallexp::Float64             = 1.5
end

# ---------------------------------------------------------------------------
# Axons
# ---------------------------------------------------------------------------

"""
    AxonParams(; kwargs...)

Parameters for axon-process generation. Defaults follow
`check_axon_params.m`.
"""
Base.@kwdef mutable struct AxonParams
    flag::Int            = 1
    distsc::Float64      = 0.5
    fillweight::Float64  = 100.0
    maxlength::Float64   = 200.0
    minlength::Float64   = 10.0
    maxdist::Float64     = 100.0
    maxel::Int           = 8
    varfill::Float64     = 0.3
    maxvoxel::Int        = 6
    padsize::Int         = 20
    numbranches::Int     = 20
    varbranches::Float64 = 5.0
    maxfill::Float64     = 0.5
    N_proc::Int          = 10
    l::Float64           = 25.0
    rho::Float64         = 0.1
end

# ---------------------------------------------------------------------------
# Background (neuropil) dendrites
# ---------------------------------------------------------------------------

"""
    BackgroundParams(; kwargs...)

Parameters for the neuropil background-dendrite generation. Defaults follow
`check_bg_params.m`.
"""
Base.@kwdef mutable struct BackgroundParams
    flag::Int           = 1
    distsc::Float64     = 0.5
    fillweight::Float64 = 100.0
    maxlength::Float64  = 200.0
    minlength::Float64  = 10.0
    maxdist::Float64    = 100.0
    maxel::Int          = 1
end

# ---------------------------------------------------------------------------
# Spike options
# ---------------------------------------------------------------------------

"""
    SpikeOptions(; kwargs...)

Top-level options for time-trace generation. Enum-like fields are stored as
`Symbol`s. Defaults follow `check_spike_opts.m`.

- `dyn_type ∈ (:AR1, :AR2, :single, :Ca_DE)` — calcium dynamics flavour.
- `prot ∈ (:GCaMP6, :GCaMP6f, :GCaMP6s, :GCaMP3, :GCaMP7)` — indicator protein.
- `rate_dist ∈ (:gamma, :uniform)` — distribution of per-neuron rates.
- `smod_flag ∈ (:hawkes, :independent)` — spike-correlation model.
"""
Base.@kwdef mutable struct SpikeOptions
    K::Int                   = 30
    mu::Float64              = 0.0
    sig::Float64             = 1.0
    dyn_type::Symbol         = :Ca_DE
    rate_dist::Symbol        = :gamma
    dt::Float64              = 1.0 / 30.0
    nt::Int                  = 1000
    rate::Float64            = 1e-3
    N_bg::Int                = 0
    prot::Symbol             = :GCaMP6
    alpha::Float64           = 1.0
    burst_mean::Float64      = 10.0
    smod_flag::Symbol        = :hawkes
    p_off::Float64           = 0.2
    selfact::Float64         = 1.2
    min_mod::Vector{Float64} = [0.4, 2.53]
    spikeflag::Bool          = true
    dendflag::Bool           = true
    axonflag::Bool           = true
end

# ---------------------------------------------------------------------------
# Calcium dynamics
# ---------------------------------------------------------------------------

"""
    CalciumParams(prot::Symbol = :GCaMP6; kwargs...)

Calcium-dynamics parameters. The four kinetic constants (`ca_amp`, `t_on`,
`t_off`, `ext_rate`) have protein-specific defaults drawn from upstream
`check_cal_params.m`; any of them can be overridden by passing the matching
keyword. `prot` is one of `:GCaMP6` (== `:GCaMP6f`), `:GCaMP6s`, `:GCaMP3`,
`:GCaMP7`; unknown symbols fall back to the `:GCaMP6` defaults.

Concentrations are in mol/L (`ca_dis`, `ca_rest`, `ind_con`) inherited
verbatim from upstream.
"""
Base.@kwdef mutable struct CalciumParams
    ca_bind::Float64  = 110.0
    ca_rest::Float64  = 50e-9
    ind_con::Float64  = 200e-6
    ca_dis::Float64   = 290e-9
    ca_sat::Float64   = 1.0
    sat_type::Symbol  = :double
    dt::Float64       = 1.0 / 30.0
    a_bind::Float64   = 3.5
    a_ubind::Float64  = 7.0
    ca_amp::Float64   = 76.1251         # :GCaMP6 / :GCaMP6f default
    t_on::Float64     = 0.8535          # :GCaMP6 / :GCaMP6f default
    t_off::Float64    = 98.6173         # :GCaMP6 / :GCaMP6f default
    ext_rate::Float64 = 292.3           # :GCaMP6 / :GCaMP6f default
    prot::Symbol      = :GCaMP6
end

function CalciumParams(prot::Symbol; kwargs...)
    cp = CalciumParams(; prot, kwargs...)
    _apply_protein_defaults!(cp, kwargs)
    return cp
end

function _apply_protein_defaults!(cp::CalciumParams, kwargs)
    p = cp.prot
    supplied = keys(kwargs)
    if p === :GCaMP6 || p === :GCaMP6f
        :ca_amp   in supplied || (cp.ca_amp   = 76.1251)
        :t_on     in supplied || (cp.t_on     = 0.8535)
        :t_off    in supplied || (cp.t_off    = 98.6173)
        :ext_rate in supplied || (cp.ext_rate = 292.3)
    elseif p === :GCaMP6s
        :ca_amp   in supplied || (cp.ca_amp   = 54.6943)
        :t_on     in supplied || (cp.t_on     = 0.4526)
        :t_off    in supplied || (cp.t_off    = 68.5461)
        :ext_rate in supplied || (cp.ext_rate = 299.0833)
    elseif p === :GCaMP7
        :ca_amp   in supplied || (cp.ca_amp   = 230.917)
        :t_on     in supplied || (cp.t_on     = 0.020137)
        :t_off    in supplied || (cp.t_off    = 3.1295)
        :ext_rate in supplied || (cp.ext_rate = 265.73)
    elseif p === :GCaMP3
        :ca_amp   in supplied || (cp.ca_amp   = 0.05)
        :t_on     in supplied || (cp.t_on     = 1.0)
        :t_off    in supplied || (cp.t_off    = 1.0)
        :ext_rate in supplied || (cp.ext_rate = 265.73)
    end
    return cp
end

# ---------------------------------------------------------------------------
# Noise model
# ---------------------------------------------------------------------------

"""
    NoiseParams(; kwargs...)

Poisson-Gauss measurement-noise parameters. Defaults follow
`check_noise_params.m`. (Upstream has a duplicate-block bug for `bleedp` /
`bleedw` whose second branch is unreachable; the values below match what the
upstream code actually applies, namely `0.3` and `0.4`.)
"""
Base.@kwdef mutable struct NoiseParams
    mu::Float64        = 100.0
    mu0::Float64       = 0.0
    sigma::Float64     = 2300.0
    sigma0::Float64    = 2.7
    darkcount::Float64 = 0.05
    sigscale::Float64  = 2e-7
    bleedp::Float64    = 0.3
    bleedw::Float64    = 0.4
end

# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------

"""
    ScanParams(; kwargs...)

Scanning parameters. Defaults follow `check_scan_params.m`.
"""
Base.@kwdef mutable struct ScanParams
    scan_buff::Int = 10
    motion::Bool   = true
    scan_avg::Int  = 2
    sfrac::Int     = 2
    verbose::Int   = 1
end

# ---------------------------------------------------------------------------
# Point-spread function
# ---------------------------------------------------------------------------

"""
    PSFFastMask(; kwargs...)

Sub-struct for the upstream `fastmask=true` code path; values match the
upstream `psf_params.FM` defaults.
"""
Base.@kwdef mutable struct PSFFastMask
    sampling::Float64 = 10.0
    fineSamp::Float64 = 2.0
    ss::Float64       = 1.0
end

"""
    PSFParams(; kwargs...)

Point-spread-function parameters. Defaults follow `check_psf_params.m`. The
`type` symbol selects the PSF flavour (`:gaussian` is implemented now;
`:vtwins` and `:bessel` are deferred). `scaling` symbol selects the photon-
scaling regime (`:two_photon`, `:three_photon`, `:temporal_focusing`).
"""
Base.@kwdef mutable struct PSFParams
    NA::Float64                 = 0.6
    objNA::Float64              = 0.8
    n::Float64                  = 1.35
    n_diff::Float64             = 0.02
    lambda::Float64             = 0.92
    obj_fl::Float64             = 4.5
    ss::Float64                 = 2.0
    sampling::Float64           = 50.0
    psf_sz::Vector{Float64}     = [20.0, 20.0, 50.0]
    prop_sz::Float64            = 10.0
    blur::Float64               = 3.0
    scatter_sz::Vector{Float64} = [0.51, 1.56, 4.52, 14.78]
    scatter_wt::Vector{Float64} = [0.57, 0.29, 0.19, 0.15]
    zernikeWt::Vector{Float64}  = [0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.12]
    taillength::Float64         = 50.0
    type::Symbol                = :gaussian
    scaling::Symbol             = :two_photon
    hemoabs::Float64            = 0.00674 * log(10)
    propcrop::Bool              = true
    fastmask::Bool              = true
    FM::PSFFastMask             = PSFFastMask()
end

# ---------------------------------------------------------------------------
# TPM photon-flux parameters
# ---------------------------------------------------------------------------

"""
    TPMParams(; kwargs...)

Two-photon-microscope photon-flux parameters. `phi` defaults to `NaN` and is
filled by [`finalize!`](@ref) from `nac` and `nidx` (matching the upstream
`isempty(tpm_params.phi)` branch). Defaults follow `check_tpm_params.m`.
"""
Base.@kwdef mutable struct TPMParams
    nidx::Float64   = 1.33
    nac::Float64    = 0.8
    phi::Float64    = NaN          # NaN → derive from nac/nidx in finalize!
    eta::Float64    = 0.6
    conc::Float64   = 10.0
    delta::Float64  = 2.0
    gp::Float64     = 0.588
    f::Float64      = 80.0
    tau::Float64    = 150.0
    pavg::Float64   = 40.0
    lambda::Float64 = 0.92
end

# ---------------------------------------------------------------------------
# Derived-field resolution
# ---------------------------------------------------------------------------

"""
    finalize!(p)

Fill any fields whose default is a sentinel (`0`, `NaN`) with values derived
from the other fields, in place. Returns `p` for chaining.

- `VolumeParams`: resolves `N_neur` from `neur_density × prod(vol_sz)` and
  `N_den` from `AD_density × xy-area`. Also rounds `vol_sz[3]` up to the
  next multiple of 10 (upstream invariant for dendrite indexing).
- `TPMParams`: resolves `phi` from `nac` / `nidx` when `isnan(phi)`.

Other parameter structs are returned unchanged.
"""
function finalize! end

function finalize!(vp::VolumeParams)
    length(vp.vol_sz) == 3 || throw(ArgumentError("vol_sz must have 3 elements"))
    if mod(vp.vol_sz[3], 10) != 0
        vp.vol_sz[3] = 10 * cld(vp.vol_sz[3], 10)
    end
    vol_um3 = prod(vp.vol_sz)
    xy_mm2  = prod(@view vp.vol_sz[1:2]) / 1e6
    if vp.N_neur == 0
        vp.N_neur = ceil(Int, vp.neur_density * vol_um3 / 1e9)
    else
        vp.neur_density = 1e9 * vp.N_neur / vol_um3
    end
    if vp.N_den == 0.0
        vp.N_den = vp.AD_density * xy_mm2
    end
    return vp
end

function finalize!(tp::TPMParams)
    if isnan(tp.phi)
        ratio = tp.nac / tp.nidx
        ratio < 1 || throw(ArgumentError("nac/nidx must be < 1 to derive phi"))
        tp.phi = 0.8 * ((1 - sqrt(1 - ratio^2)) / 2) * 0.4
    end
    return tp
end

# Identity fallback for structs without derived fields.
finalize!(p::NeuronParams)            = p
finalize!(p::VasculatureParams)       = p
finalize!(p::VasculatureNodeParams)   = p
finalize!(p::DendriteParams)          = p
finalize!(p::AxonParams)              = p
finalize!(p::BackgroundParams)        = p
finalize!(p::SpikeOptions)            = p
finalize!(p::CalciumParams)           = p
finalize!(p::NoiseParams)             = p
finalize!(p::ScanParams)              = p
finalize!(p::PSFParams)               = p
finalize!(p::PSFFastMask)             = p
