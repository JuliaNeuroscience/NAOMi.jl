# Session Handoff — 2026-05-15

## Project maturity target

`releasable-package` — `NAOMi` (UUID `04116bad-86ec-414e-852a-2781886c9025`),
targeted for the JuliaNeuroscience GitHub org.

## What was just completed

**Chunk 15 — Scanning III: full scan + motion.** Ported into
`src/scanning/scan.jl`:

- `scan_volume` — full multi-frame scan loop with AR-1-like motion
  jitter and Poisson-Gauss noise application.
- `img_sub_row_shift` — per-row fractional x/y shift + buffer crop.

End-to-end (volume → 3-frame movie with noise) runs in <1 s on the
30×30×20 µm test volume.

633 tests pass on Julia 1.10 LTS; +11 over Chunk 14.

## Key decisions made

- **`blurredBackComp2.m` not ported.** Depends on `psfT`/`psfB`
  temporal-focusing scattering background, which requires the
  cortical-light-path orchestrator deferred from Chunk 7.
- **TIFF streaming output paths skipped.** Chunk 17 will add them.
- **Motion model preserves upstream's small-step sampler.** Per-frame
  Bernoulli jumps + small uniform jitter + per-row shear vector; the
  3×nt history returns via `return_motion=true`.

## State of the codebase

- Files created or modified:
  - `src/scanning/scan.jl` — populated (was placeholder; +210 LOC).
  - `test/scanning/test_scan.jl` — new (+11 tests).
  - `test/runtests.jl` — includes the new test file.
  - `ANALYSIS_PLAN.md` — chunk-status table updated + chunk-15
    notes/deviations + ledger entry.
- Package loads cleanly: yes.
- Test suite passes: yes — 633/633 on Julia 1.10 LTS.
- Entry point(s): **the full scanning pipeline now works end-to-end**.
  A user can call `simulate_neural_volume` → `scan_volume` and obtain
  a movie.

## Next chunk

**Chunk 16 — Ideal components + ground truth.** Port
`calculateIdealComps.m`, `scan_ideal.m`, `times_from_profs.m`,
`comps2ideals.m`, `constrainEstToSomas.m`. Target file:
`src/scanning/ideal.jl`. These produce the "ideal" per-cell extracted
fluorescence traces from a clean movie + soma masks — used to compare
against downstream analysis pipelines (CNMF, Suite2P, etc).

Tests should cover: ideal traces extracted from a synthetic movie and
ideal profiles recover `neur_act.soma` to high correlation.

## Watch out for

- **`scan_volume`'s clean output** is the natural input to the ideal
  components chunk. Use `return_clean=true` to obtain `mov_clean`.
- **`set_cell_fluorescence`'s `gp_vals[k].is_soma` bitvector** marks
  which voxels are soma vs dendrite — the ideal-profile code needs
  this distinction.
- **Per-cell soma mask = projection of soma voxels onto the imaging
  plane.** The chunk-16 code should pre-compute these 2-D masks from
  the 3-D soma indices.

## Working stance reminder

`ANALYSIS_PLAN.md` "Working stance" authorises autonomous chunk
progression with auto-commits.

## Suggested next workflow

`/clear` and re-run `/new-analysis-implement`.
