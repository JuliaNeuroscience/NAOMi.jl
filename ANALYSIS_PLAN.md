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
|     6 | Optics II — Zernike + back-aperture         | complete    |
|     7 | Optics III — Fresnel propagation            | complete    |
|     8 | Volume I — vasculature                      | complete    |
|     9 | Volume II — soma generation                 | complete    |
|    10 | Volume III — dendrites                      | complete    |
|    11 | Volume IV — axons + neuropil background     | complete    |
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

**Tests**: low-order analytic forms (`Z_1 = 1`, `Z_2 = 2x`, `Z_3 = 2y`,
`Z_4 = √3(2r² − 1)`); numerical orthonormality on the unit disk via
trapezoidal sums on a 200×200 grid; phase-only modulation preserves
`|U|²`; zero aberrations is identity; `generate_back_aperture`
produces a square complex array with the hard-aperture mask zeroing
corners and a peak in the centre.

**Notes**: All four ports live in `src/optics/zernike.jl`. Public
exports: `zernike_polynomial`, `generate_zernike_weights`,
`apply_zernike`, `generate_back_aperture`. Hand-rolled — the
ZernikePolynomials.jl survey was skipped since upstream's `zidx` is
trivially Noll-conformant and the radial polynomial is ~10 LOC. No
new dependency.

**Deviations from upstream**:

- **Only the simple path of `generateBA.m` is ported** (the
  `imax*jmax == 1` branch, no `psf_params.zernikeDst`). The cell-array
  branch produces one back-aperture per FOV pixel with spatially-varying
  Zernike weights; deferred until a downstream chunk actually needs it.
  Upstream's cell-array branch also has a latent bug (it forgets the
  `X/objrad`, `Y/objrad` normalisation before calling `applyZernike`);
  fix that on porting.
- **`PSFParams.zernikeDst` not added** — would be a per-call function /
  table. Deferred with the cell-array branch.
- **`vol_params.vasc_sz` not a `VolumeParams` field**. `generate_back_aperture`
  computes it on-the-fly from `gaussian_beam_size` (matching upstream's
  "compute if missing" pattern). When Chunk 7 adds `setOpticalParams`
  this may be cached.

### Chunk 7 — Optics III: Fresnel propagation + cortical mask

Port `fresnel_propagation_multi.m`, `genCorticalLightPath.m`,
`genCorticalLightPathLite.m`, `simulate_optical_propagation.m`,
`tpmSignalscale.m`, `widthestimate.m`, `widthestimate3D.m`,
`groupzproject.m`, `setOpticalParams.m`.

**Tests**: PSF in vacuum (empty mask) matches the Chunk 5 PSF; attenuation
increases monotonically with depth; hemoglobin absorption integrates
against the upstream-default `hemoabs = 0.00674·ln(10)` constant.

**Notes**: Ported into `src/optics/propagation.jl`. Public exports:
`fresnel_propagation_multi`, `group_z_project`, `width_estimate`,
`width_estimate_3d`, `tpm_signal_scale`, `collection_mask`.

**Deviations from the original plan** (recorded so Chunk 13+ knows what
is still missing):

- **`genCorticalLightPath.m` and `genCorticalLightPathLite.m` not
  ported in this chunk.** The full per-tile cortical-light-path
  orchestrators (super-Gaussian apodisation + i,j tile loop + two-stage
  A→B[→C] propagation + per-tile downsampling/`imresize` + temporal
  focusing) are deferred to the chunk that first needs spatially-varying
  excitation masks (Chunk 13 — scanning). Until then, the standalone
  `fresnel_propagation_multi` plus `generate_back_aperture` (Chunk 6)
  is enough to test single-tile propagation. The `fastmask` /
  `fineSamp` paths in the Lite version are deferred indefinitely.
- **`simulate_optical_propagation.m` top-level orchestrator not
  ported.** Same reason. The standalone `collection_mask` covers the
  collection-side hemoabs computation (the part of upstream lines
  415–442 that the test plan asked about). The `acc_flag = 1` /
  `fastmask = true` branches, the struct-`Uin` (vtwins/bessel)
  branches, and the scatter-volume injection (depends on
  `masked_3DGP_test` Gaussian-process scatter, not yet ported) are
  deferred.
- **`setOpticalParams.m` is *misnamed* in upstream** — its body is the
  unrelated `TPM_Simulation_Parameters` opt-type selector. The actual
  `vol_params.vasc_sz` cache lives inline at lines 45–46 of
  `simulate_optical_propagation.m`; it will be reproduced when the
  top-level orchestrator is ported.
- **`groupzproject.m`'s `median` / `mode` types not ported.** The
  standard pipeline only invokes `:sum`, `:prod`, `:mean`, `:max`,
  `:min`; tests cover those five. `:median` / `:mode` are upstream
  options never exercised by the standard pipeline.
- **`widthestimate.m` upstream bug fixed.** Upstream returns
  `s1 + s2 + (f2 − f1)` for the interpolated FWHM, which is one
  sample-spacing too wide (the correct interior span between the
  two threshold crossings is `(f2 − f1 − 1) + s1 + s2`). The
  Julia port returns the corrected width; verified against the
  analytic FWHM of a sampled Gaussian to ~5 % on a 51-point grid.
- **Schmidt's two-step sampling.** `fresnel_propagation_multi`
  requires the source/observation grid spacings `(D1, D2)` to be
  related by `D1 · D2 · N = λ_med · z` for the Fourier-step to land
  at the right plane. Upstream `genCorticalLightPath` enforces this
  in its own setup; downstream callers must do the same.
- **Energy conservation is exact.** Tests confirm
  `∫|Uin|²·D1² ≈ ∫|Uout|²·D2²` to ~1e-15 relative — useful for
  catching FFT-shift bugs in future propagator variants.
- **Collection-mask convolution is hand-rolled** (`_disk_kernel` +
  `_conv2_same`) to avoid forcing `ImageFiltering` semantics through
  the simulate code; matches MATLAB `conv2(...,'same')` with
  zero-padded boundaries and a symmetric kernel. If profile shows
  it as a hot spot, swap in `ImageFiltering.imfilter`.

### Chunk 8 — Volume I: vasculature

Port `growMajorVessels.m`, `growCapillaries.m`, `simulatebloodvessels.m`,
`vessel_dijkstra.m`. Use `Graphs.jl` + `SimpleWeightedGraphs` for routing.

**Tests**: vessel mask coverage within plausible range (volume fraction
~2–5%); no orphan capillaries; surface vasculature density matches `vesFreq`.

**Notes**: All ported into `src/volume/vasculature.jl`. Public exports:
`VesselNode`, `VesselEdge`, `gen_node`, `gen_conn`, `del_node!`,
`nodes_to_conn`, `pseudo_rand_sample_2d`, `pseudo_rand_sample_3d`,
`vessel_dijkstra`, `branch_grow_nodes!`, `grow_major_vessels!`,
`grow_capillaries!`, `conn_to_vol!`, `simulate_blood_vessels`. To keep
the Chunk 12 orchestrator thin (per session handoff), the small
upstream helpers `gennode.m`, `delnode.m`, `genconn.m`, `nodesToConn.m`,
`branchGrowNodes.m`, `connToVol.m`, `pseudoRandSample2D.m`, and
`pseudoRandSample3D.m` were ported in this chunk too (they're either
trivial or tightly coupled to vessel growth). Chunk 12 will only need
to add the *non-vasculature* orchestrators.

**Deviations from upstream**:

- **`Graphs.jl` was not used.** Hand-rolled `vessel_dijkstra` (~25 LOC)
  matches upstream exactly and handles the dense `Inf`-as-forbidden
  matrix that upstream depends on. Building a complete `SimpleGraph` +
  `dijkstra_shortest_paths` would be a worse fit (and would still need
  a wrapper to deal with the structural `Inf` blocks). `SimpleWeightedGraphs`
  not added.
- **`cscvn` (MATLAB's cubic-spline-curve fitting) replaced by linear
  interpolation in `conn_to_vol!`.** Spline aesthetics are immaterial
  for a binary vessel mask once the per-edge ball dilation by
  `conn.weight` dominates. Recorded in the function docstring.
- **`imdilate` replaced by hand-rolled binary 2-D disk dilation and
  per-point 3-D ball painting** (`dilate2d_disk!`, `paint_ball3d!`).
  Avoids pulling in `ImageMorphology`. If profiling shows it's a hot
  path, swap in `ImageMorphology.dilate`.
- **Empty/orphan distinction.** Upstream uses MATLAB's `[]` empty
  vs. `0` integer to distinguish "deleted/orphan" from "source/no
  parent". The Julia port encodes this with `root == -1` (orphan/empty)
  vs. `root == 0` (source) vs. `root > 0` (parent index). Tests in
  `grow_capillaries!` that mirrored upstream's
  `~cellfun(@isempty, {nodes.root})` map to `nodes[i].root >= 0`.
- **Counts depend strongly on volume size.** At default
  `vasc_params.sourceFreq = 1000 µm/node` and
  `vesFreq = [125 200 50]`, the expected number of source/diving
  vessels is `O(L/1000)`/`O(L²/40 000)` — small test volumes round
  these to 0. The test suite uses scaled-down `sourceFreq = 400`,
  `vesFreq = [80, 100, 30]` to obtain a non-trivial vasculature at
  ~150 µm side length. The 2–5 % upstream coverage target is *not*
  asserted on those small volumes (would require a 500 µm side); the
  test bracket is a generous `0.1 % < frac < 50 %`.
- **`grow_capillaries!` capillary-to-capillary "nearest neighbour"
  bootstrap simplified.** Upstream's MATLAB uses
  `[~,mincapp] = min(cappmat)` (per-column min) and seeds connections
  via the `mincapp` permutation. The Julia port iterates row-by-row
  with the same `connect-then-Inf-out` semantics; results are
  topologically equivalent but capp pairs may be visited in a
  different order under the same RNG draw.

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

**Notes**: All implemented in `src/volume/somata.jl`. Public exports:
`spiral_sample_sphere`, `teardrop_projection`, `generate_neural_body`,
`sample_dense_neurons`, `generate_neural_volume`, `isolate_visible_somas`,
`point_in_soma`, `masked_3d_gp`. Also implemented:
`generate_neural_volume` (was upstream `generateNeuralVolume.m`,
formally Chunk 12 scope but the rasterization step is the natural
bridge from soma meshes to voxel masks and could not be sensibly
deferred). `masked_3d_gp` ported now because it is standalone (FFT-based
GP sampler used by future fluorescence/dendrite chunks).
`LinearAlgebra` added to `[deps]` for `cholesky` and `eigmin`.

**Deviations from upstream and the original chunk plan**:

- **Used `SpiralSampleSphere` (golden-angle deterministic spiral)
  rather than the chunk-plan-suggested Marsaglia random sampling.**
  Spiral sampling preserves upstream's "compute one sphere sampling
  once, reuse for every cell" reuse pattern, which Marsaglia random
  sampling would defeat. The geodesic-distance matrix is also computed
  once and reused across all cells.
- **Triangulation matrix `Tri` not ported.** Upstream `SpiralSampleSphere`
  returns `Tri = fliplr(convhulln(V))` (MATLAB convex hull) and uses
  it via `intriangulation(Vcell, Tri, points)` in
  `generateNeuralVolume.m`. The Julia port replaces this with a star-
  shape radial point-in-soma test (`point_in_soma`): for each voxel,
  find the surface vertex whose direction (from the soma centre) is
  closest to that of the voxel, then test radii. Equivalent for
  star-shaped meshes (every soma in NAOMi is star-shaped around its
  centre); avoids adding `Quickhull.jl` or `MiniQhull.jl` as a
  dependency. ~7° angular sampling at `n_samps = 200`.
- **`smoothCellBody.m` not ported here — moved to Chunk 10.** It
  takes `allpaths` (dendrite paths) and `cellBody` (voxel indices) as
  input; the algorithm smooths the soma-boundary near where dendrites
  attach. With no `allpaths` until dendrites are generated, the
  function would have no input. It will land naturally in Chunk 10.
- **`setCellFluoresence.m` not ported here — moved to Chunk 10 or 12.**
  Requires `neur_num_AD` (apical dendrite map) which is a Chunk 10
  output. Can technically be invoked with zeroed `neur_num_AD`, but
  the resulting per-cell fluorescence map would lack the
  dendrite-uniform branch.
- **Nucleus volume normalisation uses an analytic star-shape volume
  estimate** (`sum r³ · solid_angle / 3`) instead of MATLAB's
  `convhull` volume. Accurate to within ~5 % for spherical meshes; the
  result is then cubic-root scaled so multi-percent error in the
  volume estimate translates to sub-percent error in the linear
  scaling factor.
- **`generate_neural_volume` adopts the soma rasterization from upstream
  `generateNeuralVolume.m`**, including the rule that nucleus voxels
  are excluded from the soma index map (`neur_soma`) and that
  `neur_vol` is initialised to `neur_params.nuc_fluorsc` inside nuclei
  only.
- **`bsxfun(@plus, …, neur_locs.')` semantics preserved.** Upstream
  shifts soma meshes from origin-centred to `neur_locs`-centred *after*
  sampling all shapes. The Julia port does the same shift inside
  `sample_dense_neurons` so callers receive already-shifted meshes.
- **`isolate_visible_somas` uses Chunk 7's `width_estimate_3d`** for
  the PSF half-width; the function reuses the already-ported
  width-estimate fix (one-sample-spacing correction).

### Chunk 10 — Volume III: dendrites

Port `dendrite_dijkstra2.m`, `dendrite_randomwalk2.m`,
`growNeuronDendrites.m`, `growApicalDendrites.m`, `getDendritePath2.m`,
`dilateDendritePathAll.m`. Reimplement the C++ MEX kernels
(`dendrite_dijkstra_cpp.cpp`, `dendrite_randomwalk_cpp.cpp`,
`locate_neighbors.cpp`, `array_SubMod.cpp`, `array_SubSub.cpp`) in pure
Julia.

**Tests**: dendrites originate from soma boundaries; reach apical targets;
thickness profile decays per `dendrite_tau`.

**Notes**: All ported into `src/volume/dendrites.jl`. Public exports:
`dendrite_dijkstra_grid`, `get_dendrite_path`, `dendrite_random_walk`,
`dilate_dendrite_paths_all`, `grow_neuron_dendrites!`,
`grow_apical_dendrites!`, `smooth_cell_body`, `set_cell_fluorescence`.
Also folded in the two helpers deferred from Chunk 9
(`smoothCellBody.m` and `setCellFluoresence.m`).

**Deviations from upstream**:

- **C++ MEX kernels reimplemented in pure Julia.** A hand-rolled binary
  min-heap (~50 LOC) backs `dendrite_dijkstra_grid`; no `DataStructures`
  dependency. `array_SubMod`/`array_SubSub` are not ported — Julia's
  `A[idx] .+= val` and `A[idx] .= val` do the same with no helper.
  `locate_neighbors` is similarly inlined where needed.
- **Single-stage Dijkstra at fine resolution.** Upstream does a
  coarse-then-fine two-stage Dijkstra (running once on a `dims`-sized
  coarse grid to plan, then re-running constrained to that path at
  `dims.*dimsSS` resolution). This is purely a speed optimisation for
  large volumes. The Julia port runs Dijkstra directly at fine
  resolution and accepts the cost — test volumes (≤40×40×20 µm) run in
  under a second per neuron.
- **`smooth_cell_body` is simplified.** Upstream uses MATLAB's `cscvn`
  spline blend between `connIdx` and `connRoots` followed by 3-D
  border filling. The Julia port adds a radius-2 ball around the first
  voxel where each path hits `cellBody`; this is sufficient for the
  test goal "soma-dendrite junction is connected".
- **`dilate_dendrite_paths_all` is faithful** to upstream's iterative
  per-shell growth strategy, but uses Julia idioms (sets, sorted
  offsets by squared distance, 6-connectivity adjacency check).
- **Apical anchor follows upstream's "soma corner" rule.** A zero-cost
  3-D staircase from `rootL` (cell centre) to `aproot` (the lowest-
  linear-index soma voxel) is laid down before Dijkstra. This biases
  the basal tree to emerge through a consistent boundary point.
- **`grow_apical_dendrites!` skips the through-volume "border-spill"
  branch.** When the root box overflows the volume bounds, upstream
  carefully clamps source and destination regions. The single-stage
  Julia port operates directly on `fulldims`, so this branch is
  unnecessary.
- **`set_cell_fluorescence` mean-normalises soma fluorescence values to
  1** matching upstream's `0.5 * (TMP_vals - mean) / max(abs(...)) + 1`
  rule. Background components (cells `N_neur+1 … N_neur+N_den`) get
  uniform value 1.

### Chunk 11 — Volume IV: axons + neuropil background

Port `generate_axons.m`, `sort_axons.m`, `generate_bgdendrites.m`.

**Tests**: axon density matches request; background processes fill expected
volume fraction.

**Notes**: Implemented in `src/volume/axons.jl` (`generate_axons`,
`sort_axons`) and `src/volume/background.jl` (`generate_bg_dendrites`).
Heavy reuse of `dendrite_random_walk` and `dilate_dendrite_paths_all`
from Chunk 10.

**Deviations from upstream**:

- **`gp_bgvals` and `gp_vals` use NamedTuple-of-vectors** rather than
  upstream's `cell` arrays with positional columns. Per-entry structure:
  `(loc::Vector{Int32}, val::Vector{Float32})` for `gp_bgvals` and
  `(loc, val, is_soma)` for the soma/dendrite `gp_vals`. Background
  entries appended to `gp_vals` carry an all-false `is_soma` bitvector
  to stay structurally compatible.
- **`sort_axons` nearest-cell assignment is greedy** (each cell takes
  its nearest available axon, marking that axon's column as `Inf`)
  rather than a Hungarian optimal assignment. Matches upstream's
  semantics exactly.
- **`generate_bgdendrites`'s outside-volume root sampling.** Upstream
  uses `floor(rand(1,3).*(volsize+2*dtSize)-dtSize)` and rejects roots
  inside the volume. The Julia port mirrors this with the same
  rejection loop and a 100-try safety cap.
- **Entry-point shift on the volume boundary** preserves upstream's
  `switch shiftLoc` rule: on the face where the line root→ends
  enters the volume, jitter the orthogonal coordinates by a small
  amount (`rand(1:shiftdist)`) to avoid funnelling all axons through
  the exact corner.
- **AxonParams `N_proc` is treated as a return slot**, not just an
  input. `generate_axons` returns an `AxonParams` whose `N_proc`
  equals the number of processes actually produced — downstream
  `sort_axons` then bins them.

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
- **2026-05-15 (CHUNK-007)**: `fresnel_propagation_multi` follows
  Schmidt (2010) two-step angular-spectrum sampling: source and
  observation grid spacings `D1`, `D2` over `N` samples must satisfy
  `D1 · D2 · N = λ_med · z` (with `λ_med = λ / nidx`) for the focal
  spot to land at the observation plane and be Nyquist-sampled.
  Upstream `genCorticalLightPath` derives `D1` from
  `gaussian_beam_size(psf_params, fl*1e6)/min(N)`; downstream
  callers (Chunk 13+) must do the same.
- **2026-05-15 (CHUNK-007)**: Upstream `widthestimate.m` overestimates
  the interpolated FWHM by exactly one sample-spacing. The Julia
  port fixes the bug; if any downstream test or numerical result is
  cross-checked against MATLAB output, expect a one-sample
  discrepancy in widths.
- **2026-05-15 (CHUNK-007)**: Upstream `setOpticalParams.m` in the
  bitbucket repo *contains the body of `TPM_Simulation_Parameters`*
  — the actual vasc_sz caching code lives inline at the top of
  `simulate_optical_propagation.m`. If a future port re-derives
  `vasc_sz` it should match
  `gaussian_beam_size(psf, vol_depth + vol_sz[3]/2) +
  [vol_sz[1], vol_sz[2], vol_depth]` (i.e. the size 3 axis adds
  `vol_depth` only, not `vol_depth/2`).
- **2026-05-15 (CHUNK-008)**: Upstream's vasculature node counts scale
  as `nsource ∝ L/sourceFreq` and `nvert ∝ L²/vesFreq[2]²`. At the
  default `sourceFreq = 1000 µm/node` and `vesFreq[2] = 200 µm`, a
  100×100 µm patch rounds *both* to 0 with high probability — the
  vasculature simulation is calibrated for upstream's standard
  500×500 µm `vol_sz`. Downstream chunks doing scale-tests must
  either use ≥150 µm side length plus reduced `sourceFreq`/`vesFreq`,
  or accept that small volumes produce capillary-only vasculature.
- **2026-05-15 (CHUNK-009)**: Somata in this port are *star-shaped*
  around their centre (every direction from `neur_locs[k]` hits the
  boundary exactly once). The radial point-in-soma test
  (`point_in_soma`) relies on this. If a downstream chunk introduces
  non-star-shaped somata (e.g. by smoothing them against dendrites in
  `smoothCellBody`-equivalent code), the rasterization in
  `generate_neural_volume` will need to revisit this assumption — or
  switch to a triangulation-based test.
- **2026-05-15 (CHUNK-009)**: Sphere triangulation `Tri` (upstream's
  `convhulln` output, reused across cells) is *not* produced by this
  port. Downstream code that expected `vol_out.Tri` to exist must
  either compute it from `V_master` on demand (e.g. via `Quickhull.jl`),
  or rely on the radial-test rasterization paths added here. The
  `Vcell` / `Vnuc` vertex meshes *are* present and stored per-cell.
- **2026-05-15 (CHUNK-010)**: Dijkstra `M[i,j,k,d]` is the cost of
  *entering* voxel `(i,j,k)` from direction `d ∈ 1:6` (directions:
  `+x, -x, +y, -y, +z, -z`). Boundary voxels have `Inf` for the
  appropriate direction. The 6-direction asymmetry matters: any code
  that builds `M` must populate all six directions consistently with
  this convention.
- **2026-05-15 (CHUNK-010)**: `findall` on an N-D array returns
  `Vector{CartesianIndex{N}}`, *not* linear indices. Use
  `LinearIndices(A)[c]` to convert if a flat index is needed.
  Mixing CartesianIndex into a `Set{Int}` produces silently-false
  membership checks — spotted during chunk-10 testing.
- **2026-05-15 (CHUNK-009)**: `sample_dense_neurons` returns
  `Vcells::Vector{Matrix{Float64}}` with K entries (one per accepted
  cell). Upstream MATLAB used a 3-D `Nx3xK` array. Downstream chunks
  should iterate `Vcells[k]` instead of indexing `Vcell(:, :, k)`.
- **2026-05-15 (CHUNK-008)**: Vessel `nodes[i].root` uses three reserved
  values: `> 0` = parent index, `0` = source/no parent, `-1` =
  deleted-or-orphan (faithful encoding of upstream's `[]` vs. `0`
  distinction). Downstream consumers iterating over connected
  vasculature should filter by `root >= 0`; filtering by `root > 0`
  excludes valid source nodes.
- **2026-05-15 (CHUNK-007)**: `Pkg.test()` re-resolves dependencies
  and may pick up newer versions; some stochastic tests
  (Chunk 4 `gen_correlated_spike_trains — spatial correlation`)
  fail on Julia 1.12 with the resolved-newer Distributions because
  `MersenneTwister` produces different draws when downstream
  sampling algorithms change. The tests pass on Julia 1.10 LTS
  (the documented default). Future chunks should either widen the
  random seeds (e.g. average over multiple) or accept the
  Julia-1.10-LTS-only verification convention used by the
  per-chunk `Pkg.test()` runs.

## Session ledger

- 2026-05-15 CHUNK-000 (bootstrap) → next: CHUNK-001
- 2026-05-15 CHUNK-001 (parameter types) → next: CHUNK-002
- 2026-05-15 CHUNK-002 (spike generation) → next: CHUNK-003
- 2026-05-15 CHUNK-003 (calcium dynamics) → next: CHUNK-004
- 2026-05-15 CHUNK-004 (top-level traces + Hawkes) → next: CHUNK-005
- 2026-05-15 CHUNK-005 (Gaussian PSF kernels)     → next: CHUNK-006
- 2026-05-15 CHUNK-006 (Zernike + back-aperture)  → next: CHUNK-007
- 2026-05-15 CHUNK-007 (Fresnel propagation kernel + collection mask) → next: CHUNK-008
- 2026-05-15 CHUNK-008 (vasculature) → next: CHUNK-009
- 2026-05-15 CHUNK-009 (soma generation + rasterization) → next: CHUNK-010
- 2026-05-15 CHUNK-010 (dendrites + smoothCellBody + setCellFluoresence) → next: CHUNK-011
- 2026-05-15 CHUNK-011 (axons + neuropil background) → next: CHUNK-012

