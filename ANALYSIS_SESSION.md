# Session Handoff — 2026-05-15

## Project maturity target

`releasable-package` — `NAOMi` (UUID `04116bad-86ec-414e-852a-2781886c9025`),
targeted for the JuliaNeuroscience GitHub org.

## What was just completed

**Chunk 14 — Scanning II: noise model.** Ported into
`src/scanning/noise.jl`:

- `poisson_gauss_noise` — per-frame Poisson + lognormal + Gaussian
  measurement noise.
- `pixel_bleed` — per-pixel electronic bleed-through to the next
  scanned pixel.
- `apply_noise_model` — frame-by-frame movie application.

Tests verify empirical mean and variance against analytic Poisson-
Gauss predictions; both match within 5 %/10 %.

622 tests pass on Julia 1.10 LTS; +16 over Chunk 13.

## Key decisions made

- **Dynode-chain branch not ported.** Upstream's `applyNoiseModel.m`
  has a separate `noise_params.type == 'dynode'` path; the standard
  pipeline never uses it.
- **`Distributions.Poisson` used** for sampling; `exp(μ2 + σ2·randn)`
  replaces MATLAB `lognrnd`. No new deps.

## State of the codebase

- Files created or modified:
  - `src/scanning/noise.jl` — populated (was placeholder; +110 LOC).
  - `test/scanning/test_noise.jl` — new (+16 tests).
  - `test/runtests.jl` — includes the new test file.
  - `ANALYSIS_PLAN.md` — chunk-status table updated + chunk-14
    notes/deviations + ledger entry.
- Package loads cleanly: yes.
- Test suite passes: yes — 622/622 on Julia 1.10 LTS.
- Entry point(s): volume + single frame + noise all wired.

## Next chunk

**Chunk 15 — Scanning III: full scan + motion.** Port `scan_volume.m`,
`imgSubRowShift.m`, `blurredBackComp2.m`. The motion model uses AR-1
jitter per upstream. Target file: `src/scanning/scan.jl`.

Tests should cover: end-to-end produces a `(H, W, T)` movie of the
right shape; mean intensity stable across frames; motion off → frames
identical; motion on → frame-to-frame shifts bounded.

## Watch out for

- **`scan_volume_frame` is already in place** (Chunk 13). Chunk 15's
  `scan_volume!` is the outer loop that calls `scan_volume_frame` per
  frame, applies motion, and applies the noise model.
- **`imgSubRowShift` applies a per-row x-offset** — small Bresenham-
  style shifts to model fast-axis tissue motion. Hand-rolled is fine
  (~30 LOC).
- **`blurredBackComp2` is the temporal-focusing scattering background
  blur** that lives in `scan_volume_frame`. It depends on `psfT`/`psfB`
  which are not yet ported (cortical-light-path orchestrator).
  Chunk 15 can ship without it; flag as a deferred path.

## Working stance reminder

`ANALYSIS_PLAN.md` "Working stance" authorises autonomous chunk
progression with auto-commits.

## Suggested next workflow

`/clear` and re-run `/new-analysis-implement`.
