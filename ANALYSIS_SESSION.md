# Session Handoff — 2026-05-15

## Project maturity target

`releasable-package` — `NAOMi` (UUID `04116bad-86ec-414e-852a-2781886c9025`),
targeted for the JuliaNeuroscience GitHub org.

## What was just completed

**Chunk 13 — Scanning I: PSF FFT + single-frame scan.** Ported into
`src/scanning/psf_fft.jl`:

- `psf_fft` — 2-D FFT of a 3-D PSF with optional axial pre-summing.
- `single_scan` — convolve a 3-D fluorescence volume with the PSF
  (spatial or frequency domain) and return a 2-D image.
- `setup_scan_volume_frame` — build a `ScanVolume` struct with per-cell
  soma/dend/axon/nucleus voxel-and-value vectors plus the cached PSF
  FFT.
- `scan_volume_frame` — synthesize one scan frame; supports optional
  `tpm_params` for Xu-Webb signal scaling.
- `ScanVolume` — pre-processed scan struct.

606 tests pass on Julia 1.10 LTS; +16 over Chunk 12.

## Key decisions made

- **Temporal-focusing (psfT/psfB) and excitation-collection mask paths
  skipped.** These depend on the cortical-light-path orchestrator
  deferred from Chunk 7. When that lands, this code grows new accept
  branches (no breaking changes).
- **Activity-vector tolerance.** `_coerce_activity` right-pads per-cell
  vectors with zeros when they are shorter than the cell count
  tracked by the `ScanVolume`. This lets callers track soma activity
  for primary neurons only and treat background-dendrite cells as
  silent.
- **`nearest_small_prime` FFT-size optimisation skipped.** Julia's
  FFTW handles arbitrary sizes; this rounding was a MATLAB speed
  trick.

## State of the codebase

- Files created or modified:
  - `src/scanning/psf_fft.jl` — populated (was placeholder; +280 LOC).
  - `test/scanning/test_psf_fft.jl` — new (+16 tests).
  - `test/runtests.jl` — includes the new test file.
  - `ANALYSIS_PLAN.md` — chunk-status table updated + chunk-13
    notes/deviations + ledger entry.
- Package loads cleanly: yes.
- Test suite passes: yes — 606/606 on Julia 1.10 LTS (was 590).
- Entry point(s): `simulate_neural_volume` + `setup_scan_volume_frame` +
  `scan_volume_frame` work end-to-end. A user can now go from
  parameters → volume → single scan frame.
- Known issues: none introduced this chunk. Pre-existing Chunk-4
  Julia-1.12 RNG flake still pending.

## Next chunk

**Chunk 14 — Scanning II: noise model.** Port
`PoissonGaussNoiseModel.m`, `applyNoiseModel.m`, `pixel_bleed.m`.
Target file: `src/scanning/noise.jl`.

Tests should cover: empirical mean / variance of noise output match
analytic Poisson-Gauss predictions across a `μ` sweep.

## Watch out for

- **`tpm_signal_scale` from Chunk 7** is used by `scan_volume_frame` to
  rescale clean signal to physical photon counts. The noise chunk
  consumes this output as `μ` (Poisson mean) and applies the
  Poisson-Gauss model on top. Don't double-apply the TPM scale.
- **`NoiseParams.sigscale`** (Chunk 1) is upstream's per-pixel offset
  scaling computed from `tpm_signal_scale × spike_opts.dt × sfrac² /
  (vol_sz_xy × vres²)`. The simulator updates `noise_params.sigscale`
  inline; the Julia port should do the same.
- **`pixel_bleed.m` (Gaussian-blur of adjacent pixels)** is small —
  just a 2-D `Float32` convolution. Reuse Chunk 5's
  `gaussian_psf` if a 2-D variant is needed.

## Working stance reminder

`ANALYSIS_PLAN.md` "Working stance" authorises autonomous chunk
progression with auto-commits.

## Suggested next workflow

`/clear` and re-run `/new-analysis-implement`.
