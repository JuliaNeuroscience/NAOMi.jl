# Deferred-work inventory

This document inventories every part of upstream NAOMi-Sim that the
Julia port (Chunks 0–18) did **not** implement, plus algorithmic
simplifications worth revisiting. It is the Chunk-19 deliverable.

The target repository `JuliaNeuroscience/NAOMi.jl` does not yet exist on
GitHub, so this inventory is persisted as a document rather than as
GitHub issues. **When the repository is published, each `###` entry
below should be filed as a GitHub issue** — the headings, upstream-file
references, and "what a port needs" notes are written to drop straight
into issue bodies.

Upstream MATLAB/C++ source: <https://bitbucket.org/adamshch/naomi_sim>
(paths below are relative to its `code/` directory).

---

## A. Out-of-scope modules (never ported)

These were declared out of scope in `ANALYSIS_PLAN.md` from the start.

### A1. Interactive GUI

- **Upstream**: `GUI/@MovieSlider/`, `GUI/@gui/`.
- **Role**: MATLAB-figure-based interactive movie browser and parameter
  GUI.
- **Why deferred**: GUI work is a separate project; depends on a
  MATLAB-specific figure/callback model with no direct Julia analogue.
- **A port would need**: a `Makie.jl` (or `GLMakie`) viewer; best
  treated as a downstream companion package, not part of core NAOMi.

### A2. Variant simulation scripts

- **Upstream**: `TPM_Simulation_Script_{bessel,cylindrical,layer5,deep,
  shallow,sparse,misaligned,somaOnly,gcamp6s,anatomy,anatomy2,
  highActivity,lowActivity,Blood_Vessels,LowRam,small}.m`.
- **Role**: Each is a parameter-preset copy of the standard script for a
  specific imaging scenario.
- **Why deferred**: Only `TPM_Simulation_Script_standard.m` was in
  scope (ported to `examples/standard_pipeline.jl`). The variants are
  parameter presets, not new algorithms.
- **A port would need**: each variant becomes a small parameter-override
  example script; trivial once the optical paths they exercise (Bessel,
  cylindrical, temporal-focusing — see section B) are implemented.

### A3. Analysis and plotting helpers

- **Upstream**: `AnalysisAndPlotting/plotAlgorithmComparisons.m`,
  `scriptRealDataComparison.m`, `script_fitCalciumModel.m`.
- **Role**: Compare simulated output against analysis-algorithm results
  and against real data; fit calcium-model parameters.
- **Why deferred**: Downstream evaluation tooling, not part of the
  generator. `constrainEstToSomas.m` (section C8) belongs here too.
- **A port would need**: a separate analysis package; depends on having
  CNMF/Suite2p-style outputs to compare against.

### A4. Experimental utilities

- **Upstream**: `experimental/{fn_map,fn_structmerge,fn_switch,
  fn_timevector,gen_voltage_trace_wrapper,row,spk_calcium,spk_display,
  spk_gentrain}.m`.
- **Role**: Miscellaneous helpers and a voltage-trace experiment.
- **Why deferred**: Marked experimental upstream; not reached by the
  standard pipeline. Most `fn_*` helpers have idiomatic Julia
  equivalents already.

### A5. MEX self-test programs

- **Upstream**: `MEX/{array_SubModTest,array_SubSubTest}.cpp` and the
  compiled `.mexa64`/`.mexw64` artefacts.
- **Role**: Standalone correctness testers for the C++ array kernels.
- **Why deferred**: The kernels themselves were reimplemented in pure
  Julia (Chunk 10); the `*Test` variants are redundant — Julia's
  `A[idx] .+= val` is the operation, covered by the volume test suite.

### A6. Low-RAM volume variant

- **Upstream**: `VolumeCode/simulate_neural_volume_lowram.m`,
  `TPM_Simulation_Script_LowRam.m`.
- **Role**: A streaming, lower-memory volume generator and its driver.
- **Why deferred**: An optimisation variant of the in-scope
  `simulate_neural_volume.m`. The LowRam next-timepoint functions
  (section C2) belong here.
- **A port would need**: a streaming volume generator that materialises
  components lazily; only worth it if memory becomes a real constraint
  for large volumes.

---

## B. Optical paths not yet implemented

The PSF `type` / `scaling` symbols (`:gaussian | :vtwins | :bessel`,
`:two_photon | :three_photon | :temporal_focusing`) were designed to
dispatch additively — only `:gaussian` / `:two_photon` are implemented.

### B1. Temporal-focusing optics

- **Upstream**: `OpticsCode/applyTemporalFocusing.m`; the `psfT`/`psfB`
  scattering-background fields consumed in `scan_volume.m`,
  `scan_ideal.m`, `setup_scan_volume_frame.m`.
- **Role**: Models temporal-focusing microscopy: a depth-dependent
  scattering background above and below the focal plane.
- **Why deferred**: Requires the cortical-light-path orchestrator (B3);
  no part of the standard two-photon pipeline.
- **A port would need**: B3 first, then `blurredBackComp2.m` (C7) and
  the `psfT`/`psfB` accept-branches in the scanning code.

### B2. vTwINS, Bessel, and cylindrical beam PSFs

- **Upstream**: `OpticsCode/generateBesselBA.m`,
  `generateBesselProfile.m`, `generateVtwinsBA.m`,
  `generateCylindricalBA.m`; the `PSF.left`/`PSF.right` struct branches
  in the scanning code.
- **Role**: Alternative excitation-beam shapes (extended depth of field,
  dual-plane).
- **Why deferred**: The standard pipeline uses only the Gaussian beam.
- **A port would need**: per-type back-aperture generators that slot
  into `generate_back_aperture`'s dispatch; the scanning code already
  reserves the `:vtwins`/`:bessel` symbols.

### B3. Cortical-light-path orchestrator

- **Upstream**: `OpticsCode/genCorticalLightPath.m`,
  `genCorticalLightPathLite.m`, `simulate_optical_propagation.m`,
  `simulate_optical_propagation2.m`.
- **Role**: The top-level optics driver — super-Gaussian apodisation, an
  i,j tile loop, two/three-stage A→B[→C] Fresnel propagation, per-tile
  downsampling, and (in the Lite version) the `fastmask`/`fineSamp`
  approximations. Produces spatially-varying excitation masks.
- **Why deferred** (Chunk 7): The standalone `fresnel_propagation_multi`
  + `generate_back_aperture` cover single-tile propagation, which is
  enough for the standard pipeline's PSF. The full orchestrator is only
  needed for spatially-varying excitation and the scattering background.
- **A port would need**: the tile loop, the Schmidt two-step sampling
  setup (`D1·D2·N = λ_med·z`, already documented in the Working-knowledge
  section of the plan), and the scatter-volume injection (depends on
  `masked_3d_gp`, which *is* ported).
- **Note**: upstream `setOpticalParams.m` is misnamed — its body is the
  unrelated `TPM_Simulation_Parameters` opt-selector; the real
  `vasc_sz` caching lives inline in `simulate_optical_propagation.m`.

### B4. Spatially-varying Zernike aberrations

- **Upstream**: the cell-array branch of `OpticsCode/generateBA.m`
  (`imax*jmax > 1`, `psf_params.zernikeDst`).
- **Role**: One back-aperture per FOV pixel, with Zernike weights that
  vary across the field.
- **Why deferred** (Chunk 6): Only the simple single-back-aperture
  branch is exercised by the standard pipeline.
- **A port would need**: a `PSFParams.zernikeDst` field (a per-call
  function or table) and the per-pixel loop. **Upstream bug to fix on
  porting**: the cell-array branch forgets the `X/objrad`, `Y/objrad`
  normalisation before calling `applyZernike`.

---

## C. Partial-port deferrals (branches skipped inside ported files)

### C1. Calcium double-exponential kernel variants

- **Upstream**: the `'plus'` and `'min'` branches of
  `TimeTraceCode/mk_doub_exp_ker.m`.
- **Why deferred** (Chunk 3): the standard pipeline exercises only the
  default `'mult'` form.

### C2. LowRam next-timepoint calcium/spike functions

- **Upstream**: `genNextCalciumDynamics.m`, `genNextSpikeTimepoint.m`,
  `generateNextTimePoint.m`.
- **Why deferred** (Chunk 3): referenced *only* by
  `TPM_Simulation_Script_LowRam.m`; belongs with the LowRam port (A6).

### C3. Continuous-time marked point process

- **Upstream**: the `discrete_flag = false` path through
  `TimeTraceCode/markpointproc.m`.
- **Why deferred** (Chunk 4): the standard pipeline always uses the
  discrete approximation. `bin_spike_trains` *is* ported and tested for
  the continuous path if it is ever revived.

### C4. Polyphase resampling

- **Upstream**: MATLAB `resample` (polyphase) in the time-trace
  resampler and any future scan-rate conversion.
- **Why deferred** (Chunks 4, 13): replaced by linear interpolation
  (`_resample_to_user`). Adequate for the simulation rates in scope.
- **Revisit if**: spectral artefacts appear at slow `dt`.

### C5. Multi-batch time-trace generation

- **Upstream**: the `batch_sz < N_node` path of `generateTimeTraces.m`.
- **Why deferred** (Chunk 4): every upstream caller uses the default
  single batch.

### C6. `groupzproject` median / mode reductions

- **Upstream**: the `:median` / `:mode` types of
  `OpticsCode/groupzproject.m`.
- **Why deferred** (Chunk 7): the standard pipeline invokes only
  `:sum`, `:prod`, `:mean`, `:max`, `:min` (all ported and tested).

### C7. Temporal-focusing scattering background

- **Upstream**: `ScanningCode/blurredBackComp2.m`; the `psfT`/`psfB`/
  `colmask`/`mask` accept-branches of `setup_scan_volume_frame.m` and
  `scan_volume.m`.
- **Why deferred** (Chunks 13, 15): depends on B1/B3.

### C8. Fractional sub-sampling and small-prime FFT sizing

- **Upstream**: the `imresize` fractional-`sfrac` fallback and the
  `nearest_small_prime` factor-of-7 FFT-size bump in the scanning code.
- **Why deferred** (Chunk 13): `sfrac` integer binning (sum-and-
  subsample) is implemented; Julia's FFTW handles arbitrary sizes, so
  the prime-size bump was a MATLAB-only speed hack.

### C9. Dynode-chain noise model

- **Upstream**: `ScanningCode/DynodeNoiseModel.m`, reached via the
  `noise_params.type == 'dynode'` branch of `applyNoiseModel.m`.
- **Why deferred** (Chunk 14): the standard pipeline uses only the
  default Poisson–Gauss branch.

### C10. `scan_ideal` z-stack

- **Upstream**: `ScanningCode/scan_ideal.m`, which calls
  `single_scan_stack.m`.
- **Why deferred** (Chunk 16): `single_scan_stack.m` is **absent from
  the upstream repository** — `scan_ideal.m` is non-functional upstream
  in its full multi-z-offset form. Its working part (a `motion=0`
  baseline rescan) is subsumed by `calculate_ideal_comps`.
- **A port would need**: treating a multi-z-offset scan stack as new
  work, not a port — there is no upstream reference implementation.

### C11. `constrainEstToSomas`

- **Upstream**: `MiscCode/constrainEstToSomas.m`.
- **Role**: Subselects somatic components from a downstream-analysis
  `est` struct (`estactIdxs`, `estactideal`, `corrIdxs`, …).
- **Why deferred** (Chunk 16): operates on an `est` struct produced by
  component-matching analysis code that this port has not built; belongs
  with A3.

### C12. L1-penalised time-trace extraction

- **Upstream**: the `lambda > 0` (TFOCS `solver_L1RLS`) and
  `lambda < 0` (`linsolve`) branches of `MiscCode/times_from_profs.m`.
- **Why deferred** (Chunk 16): the Julia port implements the default
  `lambda = 0` NNLS path (projected gradient) and an unconstrained-LS
  option; `lambda > 0` throws a clear "not ported" error.
- **A port would need**: a FISTA-style proximal solver or
  `ProximalOperators.jl`.

### C13. Streaming TIFF append

- **Upstream**: `MiscCode/tifinitialize.m`, `tifappend.m`.
- **Why deferred** (Chunk 17): frame-by-frame append across separate
  calls is awkward with `TiffImages.jl` (the file handle must stay open
  across the scan loop). `write_tiff_blocks` covers the practical
  multi-file output need.

### C14. AVI movie export

- **Upstream**: `MiscCode/make_avi.m`.
- **Why deferred** (Chunk 17): renders frames through a MATLAB figure
  (`imagesc` + `getframe` + `VideoWriter`); needs a plotting backend.
  TIFF is the portable interchange format.

### C15. `.mat` workspace splitting

- **Upstream**: `MiscCode/saveSimulationParts.m`.
- **Why deferred** (Chunk 17): splits a monolithic MATLAB `.mat`
  workspace dump into part files; the Julia port never builds such a
  `.mat`, so the concept does not carry over.

### C16. `.fits` / `.mat` movie output

- **Upstream**: the `fitswrite` and `save(...,'-v7.3')` branches of
  `MiscCode/write_TPM_movie.m`.
- **Why deferred** (Chunk 17): `write_tpm_movie` supports `.tif` only;
  the other branches throw a clear "not ported" error. Would need
  `FITSIO.jl` / `MAT.jl`.

---

## D. Algorithmic simplifications (revisit if artefacts appear)

These produce correct results for the standard pipeline but diverge from
upstream's exact implementation; listed so a future maintainer knows
where fidelity was traded for simplicity or fewer dependencies.

### D1. Single-stage Dijkstra for dendrites

- **Upstream**: `VolumeCode/dendrite_dijkstra2.m` does a coarse-then-fine
  two-stage Dijkstra.
- **Port** (Chunk 10): runs Dijkstra directly at fine resolution. A
  speed simplification only; revisit for very large volumes.

### D2. Linear interpolation instead of `cscvn` splines

- **Upstream**: MATLAB `cscvn` cubic-spline curve fitting in
  `conn_to_vol!` (vasculature) and `smooth_cell_body` (somata).
- **Port** (Chunks 8, 10): linear interpolation. Spline aesthetics are
  immaterial for a binary mask after ball dilation.

### D3. Hand-rolled morphology

- **Upstream**: `imdilate` / `conv2` via the Image Processing Toolbox.
- **Port** (Chunks 7, 8): hand-rolled disk dilation, ball painting, and
  `conv2`-`same`. Avoids an `ImageMorphology` dependency; swap in
  `ImageMorphology.dilate` / `ImageFiltering.imfilter` if profiling
  shows a hot path.

### D4. Radial point-in-soma test instead of triangulation

- **Upstream**: `intriangulation(Vcell, Tri, points)` with `Tri` from
  MATLAB `convhulln`.
- **Port** (Chunk 9): a star-shape radial test (`point_in_soma`). Exact
  for star-shaped meshes (every NAOMi soma is star-shaped about its
  centre). The `Tri` triangulation is not produced; downstream code
  expecting `vol_out.Tri` must compute it on demand (e.g. via
  `Quickhull.jl`) or use the radial paths.

### D5. Greedy nearest-neighbour assignment

- **Upstream**: `sort_axons.m` cell-to-axon assignment.
- **Port** (Chunk 11): greedy (each cell takes its nearest available
  axon). Matches upstream semantics exactly — upstream is also greedy,
  not a Hungarian optimal assignment; recorded here only to flag that an
  optimal assignment was *considered and rejected as unfaithful*.

---

## Validation conventions carried by the port

- **No MATLAB-anchored fixtures.** Tests verify shape, parameter
  sweeps, and statistical/structural properties — not bit-exact
  agreement with upstream MATLAB output.
- **Julia 1.10 LTS is the verification baseline.** Some stochastic
  tests draw different `MersenneTwister` samples on Julia 1.12 when
  resolved-newer `Distributions` changes its sampling algorithms; the
  per-chunk `Pkg.test()` runs target 1.10 LTS.
