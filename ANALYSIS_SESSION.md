# Session Handoff — 2026-05-15

## Project maturity target

`releasable-package` — `NAOMi` (UUID `04116bad-86ec-414e-852a-2781886c9025`),
targeted for the JuliaNeuroscience GitHub org.

## What was just completed

**Chunk 19 — Deferred-work inventory.** Final chunk. The target
repository `JuliaNeuroscience/NAOMi.jl` does not yet exist on GitHub, so
GitHub issues could not be filed; instead the inventory was persisted as
`DEFERRED_WORK.md` at the repo root, with each `###` entry written to
drop straight into a GitHub issue body once the repo is published.

`DEFERRED_WORK.md` covers ~30 items in four categories:
- **A** — out-of-scope modules never ported (GUI, 16 variant scripts,
  analysis/plotting, experimental utilities, MEX self-tests, LowRam).
- **B** — optical paths not implemented (temporal focusing, vTwINS/
  Bessel/cylindrical PSFs, cortical-light-path orchestrator,
  spatially-varying Zernike).
- **C** — 16 branch-level deferrals harvested from every chunk's
  "Deviations from upstream" notes.
- **D** — 5 algorithmic simplifications flagged for future fidelity
  work.

No code or tests in this chunk.

## Project status: ALL CHUNKS COMPLETE

Chunks 0–19 are all `complete`. The Julia port of NAOMi-Sim's core
five-stage pipeline (TimeTraces, Optics, Volume, Scanning, I/O) is done.

- Test suite: 671/671 pass on Julia 1.10 LTS.
- Docs: `julia --project=docs docs/make.jl` builds cleanly; doctests
  pass.
- End-to-end: `julia --project examples/standard_pipeline.jl` runs the
  full volume → PSF → activity → scan → ideal-components → TIFF pipeline
  on a 30×30×20 µm volume in ~21 s.

## State of the codebase

- Files created or modified this chunk:
  - `DEFERRED_WORK.md` — new (the inventory).
  - `ANALYSIS_PLAN.md` — chunk-19 marked complete + notes + ledger;
    all 20 chunks now `complete`.
- Package loads cleanly: yes.
- Test suite passes: yes — 671/671 on Julia 1.10 LTS (unchanged; this
  chunk touched no `src/` or `test/` files).

## Suggested follow-up work

The porting plan is finished. Natural next steps, none of them part of
the plan:

1. **Publish the repository** to `JuliaNeuroscience/NAOMi.jl`, then file
   the `DEFERRED_WORK.md` entries as GitHub issues.
2. **Run the package-freshening skills** (`/freshen-package` or the
   individual `freshen-*` skills) — Aqua, ExplicitImports, struct
   mutability, coverage, gitignore — to bring the package up to
   release-quality polish.
3. **Cross-version check**: the suite is verified on Julia 1.10 LTS;
   a couple of stochastic tests are known to draw differently on 1.12
   (see the plan's Working-knowledge notes). Consider widening seeds
   or marking those tests version-tolerant before tagging a release.

## Working stance reminder

`ANALYSIS_PLAN.md` "Working stance" authorised autonomous chunk
progression with auto-commits throughout.
