# NAOMi.jl Porting Plan

This document drives `/new-analysis-implement`. Each session implements **one**
chunk, updates the "Status" column below, and prepares a handoff for the next
session.

## Status

| Chunk | Title                                       | Status      |
|------:|---------------------------------------------|-------------|
|     0 | Bootstrap                                   | complete    |
|     1 | Parameter types                             | complete    |
|     2 | TimeTraces I — spike generation             | complete    |
|     3 | TimeTraces II — calcium dynamics            | complete    |
|     4 | TimeTraces III — top-level + correlation    | complete    |
|     5 | Optics I — PSF kernels                      | complete    |
|     6 | Optics II — Zernike + back-aperture         | pending     |
|     7 | Optics III — Fresnel propagation            | pending     |
|     8 | Volume I — vasculature                      | pending     |
|     9 | Volume II — soma generation                 | pending     |
|    10 | Volume III — dendrites                      | pending     |
|    11 | Volume IV — axons + neuropil background     | pending     |
|    12 | Volume V — top-level orchestration          | pending     |
|    13 | Scanning I — PSF FFT + single-frame scan    | pending     |
|    14 | Scanning II — noise model                   | pending     |
|    15 | Scanning III — full scan + motion           | pending     |
|    16 | Ideal components + ground truth             | pending     |
|    17 | I/O + reference demo script                 | pending     |
|    18 | Documentation pass                          | pending     |
|    19 | Deferred-work inventory                     | pending     |

## Working stance

This is a port of someone else's code (the upstream MATLAB NAOMi-Sim);
the user is not in a position to easily weigh in on design decisions
mid-chunk. Default to autonomous progress:

- **Make the reasonable call and continue.** When upstream has an
  ambiguity (multiple plausible Julia idioms, an upstream bug, a name
  collision), pick the option that follows existing conventions in this
  port (preserved upstream field names, snake_case function names,
  no new dependencies unless clearly justified) and record the choice
  in the chunk's `Notes`. Don't pause to ask.
- **Auto-commit after every green chunk.** Stage chunk-relevant files
  and commit with the standard message format
  (`Chunk N: <short title>` + `Co-Authored-By: Claude Opus 4.7`).
  Never push. Never `--no-verify`.
- **Continue to the next chunk if context allows.** After committing,
  check `/context`. If free space is **> 50 %** and the just-completed
  chunk had no `blocked` issues, start the next chunk in the same
  session. Otherwise stop and let the user `/clear`.
- **Stop early on any of:** a chunk marked `blocked`; a test failure
  that isn't an obvious tolerance issue; a missing upstream file that
  needs the user to clarify scope; or context dropping below 50 % free.
- **Bias toward smaller, more frequent commits** over one big commit
  spanning multiple chunks — easier for the user to review after the
  fact.

## Context

NAOMi (Neural Anatomy and Optical Microscopy) is a simulator of two-photon
calcium-imaging volumes — anatomy + activity + optics + scanning — published by
Song, Charles, et al. (J. Neurosci. Methods, 2021). Upstream is MATLAB
(~11.4K LOC) plus C++ MEX kernels (~0.5K LOC) at
`https://bitbucket.org/adamshch/naomi_sim`.

This Julia port matters because:

- The MATLAB pipeline requires a paid toolbox chain (Image Processing, Signal,
  Statistics, plus a compiler for MEX), limiting reproducibility.
- Downstream Julia calcium-imaging packages need a Julia-native synthetic
  ground-truth source for benchmarking.
- A Julia rewrite composes cleanly with `Distributions.jl`, `Graphs.jl`,
  `ImageFiltering.jl`, `TiffImages.jl`, `FFTW.jl`.

Scope and approach:

- **Core pipeline.** Five upstream modules (Volume, Optics, TimeTraces,
  Scanning, I/O) sufficient to reproduce `TPM_Simulation_Script_standard.m`.
  Defer GUI, AnalysisAndPlotting, experimental/, MEX testers, low-RAM
  variant, and all variant scripts.
- **Statistical/functional validation only.** No MATLAB-anchored fixtures.
  Tests verify shape, parameter sweeps, statistical properties.
- **Prefer Julia equivalents** for bundled external code (S2_Sampling_Suite,
  inpaint_nans, intriangulation, sbmvivo).

## License compliance

Upstream `license.md` is MIT (Copyright 2021 Alex Song, Adam Charles). Local
`LICENSE` is MIT (Copyright 2026 Tim Holy and contributors). Compatible.
Rules for this port:

- Every Julia source file that is a direct port of an upstream `.m` file
  carries a short header citing the upstream filename and `Copyright 2021
  Alex Song, Adam Charles, MIT`.
- Do **not** vendor any upstream `.m` or `.cpp` source into this repo;
  reimplement from algorithmic descriptions and docstrings only.
- Citation (Song et al., *J. Neurosci. Methods* 358, 2021) prominent in
  README and Documenter index.
- `NOTICE.md` records upstream attribution.

## Architecture

```
src/NAOMi.jl                 # top-level module, re-exports
src/params.jl                # all @kwdef parameter structs + defaults
src/timetraces/              # spike + calcium dynamics
src/optics/                  # PSF, Zernike, Fresnel propagation
src/volume/                  # vasculature, somata, dendrites, axons, bg
src/scanning/                # scan kernel, noise, motion, ideal comps
src/io.jl                    # TIFF / AVI writers + demo script
```

Conventions:

- One file per logical group (not one-file-per-upstream-`.m`).
- Parameter structs use `@kwdef` with upstream defaults verbatim. Derived
  fields computed in a constructor or `finalize!` helper.
- Method signatures use abstract types (`AbstractArray`, `AbstractFloat`)
  per the global Julia style guide; concrete types only where layout matters.
- Functions take parameter structs by value; modifying functions use `!`
  suffix.

## Chunks

### Chunk 0 — Bootstrap

- Install this `ANALYSIS_PLAN.md` (this commit).
- Add `NOTICE.md` with upstream attribution.
- Expand `Project.toml` `[deps]` and `[compat]` for the libraries the port
  will need: `Distributions`, `StatsBase`, `Random`, `FFTW`, `AbstractFFTs`,
  `Graphs`, `ImageFiltering`, `OffsetArrays`, `TiffImages`. Bounds compatible
  with Julia 1.10.
- Run `Pkg.resolve()` via the Julia MCP server.
- Stub `src/NAOMi.jl` with submodule includes pointing at empty placeholder
  files for each future module.
- Expand `README.md` with one paragraph of scope and the upstream citation.

**Done when**: `Pkg.test()` passes (with only a placeholder test), all
submodule files exist (even if empty), Project.toml has all deps with compat
bounds, and `Pkg.resolve()` succeeds.

### Chunk 1 — Parameter types

Translate all twelve `check_*_params.m` / `check_*_opts.m` files into `@kwdef`
structs in `src/params.jl`: `VolumeParams`, `NeuronParams`, `VasculatureParams`,
`DendriteParams`, `AxonParams`, `BackgroundParams`, `SpikeOptions`,
`CalciumParams`, `NoiseParams`, `PSFParams`, `ScanParams`, `TPMParams`.
Defaults verbatim from upstream.

**Tests**: defaults match; derived fields (e.g. `N_neur` from density)
compute correctly; parameter overriding works.

**Notes**: Implemented as `Base.@kwdef mutable struct`s so derived fields
can be filled in place via `finalize!`. Two additional sub-structs created:
`VasculatureNodeParams` (the upstream `vasc_params.node_params` sub-struct)
and `PSFFastMask` (the upstream `psf_params.FM` sub-struct).
`CalciumParams(prot::Symbol; kwargs...)` constructor selects per-protein
defaults (GCaMP6/6f/6s/3/7) at construction time and lets user kwargs
override them.

**Design choices for downstream chunks**:
- **Field names preserved verbatim from upstream MATLAB** (e.g. `vesSize`,
  `dtParams`, `randWeightScale`, `objNA`, `zernikeWt`). The mixed
  camelCase/snake_case is non-idiomatic Julia, but it makes upstream
  documentation, example scripts, and Bitbucket issues map 1:1 to Julia
  code. Future ports of `simulate_*` functions can read field names
  directly from upstream `.m` source.
- **Enum-like strings → `Symbol`** (`:Ca_DE`, `:GCaMP6`, `:gaussian`,
  `:two_photon`). All `===` comparisons in downstream code should use
  Symbols, not Strings.
- **Sentinel-based derived fields**: `VolumeParams.N_neur == 0` and
  `TPMParams.phi == NaN` mean "derive me on `finalize!`". Downstream code
  that consumes these structs should call `finalize!` before reading
  derived fields; tests confirm idempotence.
- **Upstream bug noted**: `check_noise_params.m` has duplicate `bleedp` /
  `bleedw` blocks; the second is dead code. Effective upstream defaults
  (0.3, 0.4) are encoded; this is asserted in the test suite.
- **Out of scope here**: input validation beyond what upstream `check_*`
  did (e.g. range checks, symbol membership). Future chunks may add
  validation as the consuming code surfaces requirements.

### Chunk 2 — TimeTraces I: spike generation

Port `markpointproc.m`, `gen_burst_spike_times.m`, `binSpikeTrains.m`. Use
`Distributions.jl` for Gamma-distributed rates and log-normal firing
strengths.

**Tests**: empirical rate matches `SpikeOptions.rate` within tolerance;
burst-mean matches.

**Notes**: All three ported into `src/timetraces/spikes.jl`. Public
exports: `sample_firing_rates`, `generate_burst_spike_times`,
`sample_marked_point_process`, `bin_spike_trains`. `Random` is namespaced
(use `Random.default_rng()` — `default_rng` is not auto-exposed by
`using Random` in 1.10). Two-arg overloads accept an `AbstractRNG` first
for deterministic testing; no-arg forms use the default RNG.

`sample_marked_point_process` is the generic Ogata-thinning routine
needed for Chunk 4's Hawkes implementation; the dynamic-growth path
(when `nummax` is infinite) doubles the event buffer as needed. CIF
callbacks receive subarray views of the past-event arrays — do not
mutate them.

Test deps `Random` and `Statistics` added to `[extras]` /
`[targets.test]` for stochastic testing. Tests rely on fixed-seed
`MersenneTwister` for reproducibility.

### Chunk 3 — TimeTraces II: calcium dynamics

Port `make_calcium_impulse.m`, `calcium_dynamics.m`,
`mk_doub_exp_ker.m` (default branch). Hand-code GCaMP3 / GCaMP6 ODE
coefficients from the upstream constants rather than importing `sbmvivo`.

**Tests**: impulse-response decay matches expected AR poles; ODE steady
state ≈ `ca_rest` for quiet rows; impulse drives `C` above rest in every
branch; `ca_sat < 1` caps `C`; Hill curve finite/positive/monotone across
all ten supported indicator-protein symbols; unknown protein warns and
falls back to GCaMP6f.

**Notes**: All four functions ported into `src/timetraces/calcium.jl`.
Public exports: `make_doub_exp_kernel`, `make_calcium_impulse`,
`calcium_dynamics`, `fluorescence`.

**Deviations from the plan as originally drafted**:

- **`dyn_type` branches `:AR1` / `:AR2` removed from scope** — they appear
  only in the upstream *docstrings* of `generateNextTimePoint.m` and
  `generateTimeTraces.m`, never in the dispatch logic of
  `calcium_dynamics.m`. The actual `sat_type` set is
  `(:single, :Ca_DE, :double)`, all three of which are ported.
- **`genNextCalciumDynamics.m`, `genNextSpikeTimepoint.m`,
  `generateNextTimePoint.m` deferred to issue tracker (Chunk 19)** —
  `grep` confirms these are referenced *only* by
  `TPM_Simulation_Script_LowRam.m`, which is explicitly out of scope.
  Re-introducing them belongs with the rest of the LowRAM port if/when
  that ever happens.
- **`mk_doub_exp_ker.m`'s `'plus'` and `'min'` branches deferred** — the
  standard pipeline only exercises the default `'mult'` form.
- **AR-style impulse response is implemented as a hand-rolled
  difference-equation evaluator** (`_poly_from_roots` + recursion) to
  avoid a `ControlSystems.jl` dependency. Result is bit-exact against
  the analytic two-pole closed form.
- **`Ca_DE` convolution is hand-rolled** (`_conv_full_decimate`) — keeps
  `DSP.jl` off the dependency list while supporting the upstream
  `over_samp` decimation in one step.
- **`fluorescence` (upstream private `sat_nonlin`) supports all ten
  indicator-protein symbols** from `calcium_dynamics.m`. Aliases mirror
  upstream's `case {'gcamp6','gcamp6f'}`-style multi-match (e.g.
  `:OGB1`/`:OGB_1`, `:GCaMP6_RS06`/`:GCaMP6rs06`).

### Chunk 4 — TimeTraces III: top-level + correlation

Port `generateTimeTraces.m`, `genCorrelatedSpikeTrains2.m`,
`expression_variation.m`, and the supporting `sampSmallWorldMat.m`.
Hawkes correlation uses pairwise neuron distances; take locations as
an argument (the volume stage produces them; here exercised with
synthetic positions).

**Tests**: end-to-end `generate_time_traces(spike_opts)` returns the
expected `(K, nt)` shape; AR1/AR2 and Ca_DE branches both produce
finite output; Hawkes path with `N_bg > 0` produces matching shapes;
expression-variation factors are positive (or zero for silenced
cells) with the documented log-normal / uniform support; spatial
correlation: in a 1-D arrangement, pairs ≤ 2 apart have higher
50-sample-binned correlation than pairs ≥ 15 apart.

**Notes**: All four ports live in `src/timetraces/traces.jl`. Public
exports: `samp_small_world_mat`, `expression_variation`,
`gen_correlated_spike_trains`, `generate_time_traces`. The continuous
(`discrete_flag = false`) path through `markpointproc` is *not*
ported — the standard pipeline always uses the discrete approximation.

**Deviations from upstream**:

- **`extSc`/`inbSc` sized to `N_tot = K + N_bg`** rather than upstream's
  `K`. Upstream is latent-buggy when `N_bg > 0` (vector-size mismatch
  inside the discrete loop); the port fixes this without semantic
  change for `N_bg = 0`.
- **MATLAB `resample` (polyphase) replaced by linear interpolation**
  in `_resample_to_user`. Adequate for the simulation rates the
  downstream pipeline cares about; a polyphase port is left as a
  TODO if anyone observes spectral artefacts in Chunk 13+.
- **`bin_spike_trains` not invoked** by `gen_correlated_spike_trains`
  — the discrete simulation already emits at the bin grid, so the
  intermediate continuous-time events were skipped. `bin_spike_trains`
  remains exported and tested from Chunk 2 for the continuous path
  when/if that is ported.
- **Per-compartment `ext_rate` and `ext_mult` overrides preserved
  verbatim** from `generateTimeTraces.m` (`:single`/`:double` use
  hard-coded soma/dend/bg values; `:Ca_DE` overrides `bg`'s
  `ext_rate` and `dend`/`bg`'s `ext_mult`).
- **AR1/AR2 branches** *are* ported here (they live in
  `generateTimeTraces.m`, dispatching on `spike_opts.dyn_type`).
  Chunk 3's plan correctly notes that `calcium_dynamics.m` itself
  never dispatches on AR1/AR2 — but `generate_time_traces` does, via
  a separate convolution branch using `make_calcium_impulse`.
- **Multi-batch path (`batch_sz < N_node`) skipped** — every upstream
  caller uses the default single-batch.

### Chunk 5 — Optics I: PSF kernels

Port `gaussian_psf.m`, `gaussian_psf_na.m`, `gaussianBeamSize.m`,
`generateGaussianProfile.m`.

**Tests**: shape and centring of returned arrays; PSF max at origin
equals 1 (peak normalisation, *not* integral); for the NA variant, the
axial half-intensity plane sits at `z = ±psflen/2`
(where `psflen = 0.626·λ/(n − √(n² − NA²))`); `gaussian_beam_size`
returns a 2-D-only triple, monotone in `dist`; `generate_gaussian_profile`
zeroes the field outside the hard aperture.

**Notes**: All four ports live in `src/optics/psf.jl`. Public exports:
`gaussian_psf`, `gaussian_psf_na`, `gaussian_beam_size`,
`generate_gaussian_profile`. These functions are pure analytic
evaluations — no FFTs (those land in Chunk 13). A `PSFContext` struct
was suggested in the original plan; it isn't actually needed yet
(there's nothing to cache without FFT plans).

**Deviations from the original plan**:

- **PSFs are peak-normalised, not unit-integral**. Upstream
  `gaussian_psf*` returns a kernel with peak value 1.0; integration
  over the support is not unity (and is grid-dependent). Tests check
  peak normalisation instead.
- **Coordinate-axis quirk preserved verbatim**: `gaussian_psf` uses
  `(1:N) − round(N/2)` indexing (origin at `round(N/2)`) while
  `gaussian_psf_na` uses `(0:N−1) − round(N/2)` (origin at
  `round(N/2) + 1`). The two functions disagree on which sample is "zero"
  by one. Downstream code that pulls a center index must compute it
  per-function.
- **`PSFContext` deferred** until Chunk 13 wants FFT plans.

### Chunk 6 — Optics II: Zernike + back-aperture

Port `zernike.m`, `generateZernike.m`, `applyZernike.m`, `generateBA.m`.
Survey `ZernikePolynomials.jl` first; reuse if it provides Noll-indexed
evaluation, otherwise hand-roll (it's <100 LOC).

**Tests**: orthogonality on the unit disk; aberration application preserves
total power.

### Chunk 7 — Optics III: Fresnel propagation + cortical mask

Port `fresnel_propagation_multi.m`, `genCorticalLightPath.m`,
`genCorticalLightPathLite.m`, `simulate_optical_propagation.m`,
`tpmSignalscale.m`, `widthestimate.m`, `widthestimate3D.m`,
`groupzproject.m`, `setOpticalParams.m`.

**Tests**: PSF in vacuum (empty mask) matches the Chunk 5 PSF; attenuation
increases monotonically with depth; hemoglobin absorption integrates
against the upstream-default `hemoabs = 0.00674·ln(10)` constant.

### Chunk 8 — Volume I: vasculature

Port `growMajorVessels.m`, `growCapillaries.m`, `simulatebloodvessels.m`,
`vessel_dijkstra.m`. Use `Graphs.jl` + `SimpleWeightedGraphs` for routing.

**Tests**: vessel mask coverage within plausible range (volume fraction
~2–5%); no orphan capillaries; surface vasculature density matches `vesFreq`.

### Chunk 9 — Volume II: soma generation

Port `generateNeuralBody.m`, `smoothCellBody.m`, `setCellFluoresence.m`,
`pseudoRandSample2D.m`, `pseudoRandSample3D.m`, `sampleDenseNeurons.m`,
`isolateVisibleSomas.m`, `teardrop_poj.m`. Replace `S2_Sampling_Suite` with
sphere-surface sampling via uniform random rotations (Marsaglia method via
`Random.randn`). Implement GP soma-roughness directly with
Karhunen–Loève / Cholesky factorization.

**Tests**: neuron count matches `N_neur`; minimum pairwise distance
respected; eccentricity bounded by `eccen`; nucleus fluorescence ratio
matches `nuc_fluorsc`.

### Chunk 10 — Volume III: dendrites

Port `dendrite_dijkstra2.m`, `dendrite_randomwalk2.m`,
`growNeuronDendrites.m`, `growApicalDendrites.m`, `getDendritePath2.m`,
`dilateDendritePathAll.m`. Reimplement the C++ MEX kernels
(`dendrite_dijkstra_cpp.cpp`, `dendrite_randomwalk_cpp.cpp`,
`locate_neighbors.cpp`, `array_SubMod.cpp`, `array_SubSub.cpp`) in pure
Julia.

**Tests**: dendrites originate from soma boundaries; reach apical targets;
thickness profile decays per `dendrite_tau`.

### Chunk 11 — Volume IV: axons + neuropil background

Port `generate_axons.m`, `sort_axons.m`, `generate_bgdendrites.m`.

**Tests**: axon density matches request; background processes fill expected
volume fraction.

### Chunk 12 — Volume V: top-level orchestration

Port `simulate_neural_volume.m`, `branchGrowNodes.m`, `gennode.m`,
`delnode.m`, `nodesToConn.m`, `connToVol.m`, `resampVolume.m`, `genconn.m`.

**Tests**: end-to-end smoke generates a tiny 30×30×20 µm volume in <60 s
with consistent component-array dimensions.

### Chunk 13 — Scanning I: PSF FFT + single-frame scan

Port `psf_fft.m`, `single_scan.m`, `scan_volume_frame.m`,
`setup_scan_volume_frame.m`.

**Tests**: single-frame intensity at a known neuron location is positive
and increases with `TPMParams.pavg`.

### Chunk 14 — Scanning II: noise model

Port `PoissonGaussNoiseModel.m`, `applyNoiseModel.m`, `pixel_bleed.m`.

**Tests**: empirical mean/variance of noise output match analytic
Poisson-Gauss predictions across a `μ` sweep.

### Chunk 15 — Scanning III: full scan + motion

Port `scan_volume.m`, `imgSubRowShift.m`, `blurredBackComp2.m`. Motion model
uses AR-1 jitter per upstream.

**Tests**: end-to-end produces a `(H, W, T)` movie of the right shape; mean
intensity stable across frames; motion off → frames identical; motion on →
frame-to-frame shifts bounded.

### Chunk 16 — Ideal components + ground truth

Port `calculateIdealComps.m`, `scan_ideal.m`, `times_from_profs.m`,
`comps2ideals.m`, `constrainEstToSomas.m`.

**Tests**: ideal traces extracted from a synthetic movie and ideal profiles
recover `neur_act.soma` to high correlation.

### Chunk 17 — I/O + reference demo script

Port `tiff_writer.m`, `tifwrite.m`, `tifwriteblock.m`, `tifappend.m`,
`tifinitialize.m`, `tiff_reader.m`, `tifread.m`, `make_avi.m`,
`saveSimulationParts.m`, `write_TPM_movie.m`. Use `TiffImages.jl`; AVI via
`VideoIO.jl` if light, otherwise skip and document. Translate
`TPM_Simulation_Script_standard.m` to `examples/standard_pipeline.jl`.

**Tests**: round-trip TIFF write/read preserves shape and dtype; demo
script runs on a 30×30×20 µm volume in <2 min.

### Chunk 18 — Documentation pass

Docstrings for every exported symbol; Documenter.jl site under `docs/` with
sections matching the five modules; "getting started" page reproducing the
demo; citation on the index. Hook into the existing CI workflow.

### Chunk 19 — Deferred-work inventory

File GitHub issues for everything out-of-scope: GUI, all variant scripts,
low-RAM volume variant, analysis-and-plotting helpers, experimental/
utilities, MEX self-tests. Each issue references the upstream `.m` files
and the algorithmic role.

## Verification (end-to-end)

Per chunk: `Pkg.test()` from the package's environment must pass. Use the
Julia MCP server (`julia +1`) for incremental Revise-driven iteration so
compilation cost is amortized; reserve `Pkg.test()` for end-of-chunk
verification.

End-to-end (after Chunk 17): `julia --project examples/standard_pipeline.jl`
must produce a TIFF stack and a clean-TIFF stack on a 30×30×20 µm volume
within ≤2 minutes on the dev machine. Output movie must have the expected
shape, non-negative pixel values, and roughly Poisson-Gaussian intensity
statistics.

## Out of scope (issue-tracked in Chunk 19)

- GUI (`code/GUI/@MovieSlider`, `@gui`).
- Variant scripts: `_bessel`, `_cylindrical`, `_layer5`, `_deep`, `_shallow`,
  `_sparse`, `_misaligned`, `_somaOnly`, `_gcamp6s`, `_anatomy`, `_anatomy2`,
  `_highActivity`, `_lowActivity`, `_Blood_Vessels`, `_LowRam`.
- Analysis & plotting (`code/AnalysisAndPlotting/*`).
- Experimental utilities (`code/experimental/*`).
- Low-RAM volume variant (`simulate_neural_volume_lowram.m`).
- MEX self-test programs.
- Temporal focusing / vtwins / Bessel beam optics paths (the PSF type field
  carries `:gaussian | :vtwins | :bessel`; we implement only `:gaussian`
  initially and dispatch on the symbol so adding the others later is
  additive).

## Working knowledge

- **2026-05-15 (CHUNK-001)**: Parameter-struct field names are preserved
  verbatim from upstream MATLAB (mixed camelCase/snake_case) for
  traceability. Downstream chunks must not rename them in passing.
- **2026-05-15 (CHUNK-001)**: Enum-like fields are stored as `Symbol`.
  Downstream `===` dispatch should use Symbols (`:Ca_DE`, `:GCaMP6`,
  `:gaussian`, …), never Strings.
- **2026-05-15 (CHUNK-001)**: Derived fields use sentinels (`N_neur == 0`,
  `phi == NaN`) resolved by `finalize!`. Downstream consumers must call
  `finalize!` (or check for the sentinel) before reading.
- **2026-05-15 (CHUNK-001)**: Upstream `check_noise_params.m` has dead
  duplicate-block code for `bleedp` / `bleedw`; the *effective* upstream
  defaults are `0.3` / `0.4`, encoded here. If downstream tests ever
  cross-check against upstream MATLAB output, this is the value that
  matches what upstream actually applied.
- **2026-05-15 (CHUNK-002)**: `Random.default_rng` is *not* auto-exposed
  by `using Random` in Julia 1.10 — qualify it as `Random.default_rng()`
  in all submodules. (`Distributions` and other deps are imported via
  `using` in `src/NAOMi.jl` at the package top level; all submodule
  files share that namespace.)
- **2026-05-15 (CHUNK-002)**: Upstream's spike-time loop writes
  `S(k, bin) = 1` (assignment, not increment), so coincident bursts in a
  single bin do not accumulate. The Julia port mirrors this. If
  downstream code wants count semantics (e.g. for binning a Hawkes
  process), use `bin_spike_trains`, which *does* accumulate.
- **2026-05-15 (CHUNK-003)**: `CalciumParams.t_on` / `t_off` are
  *rates* (1/s), not time constants, despite the misleading names.
  This matches upstream `mk_doub_exp_ker.m`, where the kernel is
  `A·(1 − exp(−t_on·t))·exp(−t_off·t)`. With the GCaMP6f defaults
  (`t_on = 0.8535`, `t_off = 98.6173`) the kernel peaks near 50 ms and
  decays over ~150 ms. Future chunks must not "fix" this by inverting.
- **2026-05-15 (CHUNK-003)**: Quiet (`S = 0`) steady state of the
  `:single` and `:double` ODE branches is *not exactly* `ca_rest` —
  there is a small slow drift (~1 % of `ca_rest` over hundreds of
  samples) because `CB_i(0) = 0` is not the equilibrium of the binding
  ODE. This is upstream behaviour, mirrored here. Tests bound the
  drift to <1 % rather than asserting equality.
- **2026-05-15 (CHUNK-004)**: Upstream's "K_conn" parameter to
  `sampSmallWorldMat` is *not* the number of connections per node —
  the Toeplitz lattice init yields `K_conn − 1` connections per
  interior node (the diagonal is counted once), and rewiring preserves
  that count. Tests should bound row-sums against `K_conn − 1` (and
  the diagonal already carries a `1` before `self_ex` is added).
- **2026-05-15 (CHUNK-004)**: `SpikeOptions.smod_flag === :hawkes`
  generates *both* soma and background spike trains in one call
  through `gen_correlated_spike_trains`. The non-Hawkes path
  (`:independent`) generates soma and bg separately with
  `generate_burst_spike_times`. Future chunks integrating with
  generated volumes should pick the path based on `smod_flag`.
- **2026-05-15 (CHUNK-004)**: Internal simulation rate is fixed at
  100 Hz; the user-facing `spike_opts.dt` triggers a linear-interp
  resample. Polyphase / anti-aliasing is deferred; if Chunk 13+ shows
  spectral artefacts at slow `dt`, revisit.
- **2026-05-15 (CHUNK-005)**: Upstream's "axial FWHM" in
  `gaussian_psf*` is *not* the conventional Gaussian FWHM — it is the
  plane where intensity (the un-squared field amplitude) drops to ½ of
  the on-axis peak. With `psf = intensity^2`, the squared PSF at that
  plane is ¼ of its peak. `psflen` is twice this half-width.
- **2026-05-15 (CHUNK-005)**: `gaussian_psf` and `gaussian_psf_na`
  use different coordinate indexing — the former centers the origin
  on index `round(N/2)`, the latter on `round(N/2) + 1`. Downstream
  code consuming PSF arrays must derive the center per-function.

## Session ledger

- 2026-05-15 CHUNK-000 (bootstrap) → next: CHUNK-001
- 2026-05-15 CHUNK-001 (parameter types) → next: CHUNK-002
- 2026-05-15 CHUNK-002 (spike generation) → next: CHUNK-003
- 2026-05-15 CHUNK-003 (calcium dynamics) → next: CHUNK-004
- 2026-05-15 CHUNK-004 (top-level traces + Hawkes) → next: CHUNK-005
- 2026-05-15 CHUNK-005 (Gaussian PSF kernels)     → next: CHUNK-006

