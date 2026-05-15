# Session Handoff — 2026-05-15

## Project maturity target

`releasable-package` — `NAOMi` (UUID `04116bad-86ec-414e-852a-2781886c9025`),
targeted for the JuliaNeuroscience GitHub org.

## What was just completed

**Chunk 18 — Documentation pass.** Structural docs work — all 93
exported symbols already had docstrings (verified zero gaps via
`Base.Docs.doc`).

- `docs/make.jl`: site restructured into Home, Getting started, and a
  six-page API section; `checkdocs=:exports` added.
- `docs/src/index.md`: overview, five-stage table, install
  instructions, Song et al. 2021 citation with DOI.
- `docs/src/getting-started.md` (new): end-to-end walkthrough mirroring
  `examples/standard_pipeline.jl`.
- `docs/src/{parameters,timetraces,optics,volume,scanning,io}.md` (new):
  per-module API pages with explicit `@docs` blocks; every exported
  symbol listed exactly once.

The local `docs/make.jl` build succeeds and Documenter's doctest stage
runs clean. The existing `.github/workflows/CI.yml` already deploys
docs — no CI changes needed.

## Key decisions made

- **`checkdocs=:exports`.** Two internal helpers (`dilate2d_disk!`,
  `paint_ball3d!`) have docstrings but are not exported; Documenter
  1.x's default `:all` flagged them. `:exports` scopes the manual-
  inclusion check to the public API.
- **Explicit `@docs` blocks** (not `@autodocs`) so API pages control
  grouping and ordering.
- **Getting-started uses plain code blocks**, not `@example`, to keep
  the doc build fast (the real pipeline takes ~21 s).

## State of the codebase

- Files created or modified:
  - `docs/make.jl` — restructured page tree + `checkdocs=:exports`.
  - `docs/src/index.md` — rewritten (overview + citation).
  - `docs/src/getting-started.md` — new.
  - `docs/src/parameters.md`, `timetraces.md`, `optics.md`,
    `volume.md`, `scanning.md`, `io.md` — new API pages.
  - `ANALYSIS_PLAN.md` — chunk-status table + chunk-18 notes + ledger.
- Package loads cleanly: yes.
- Test suite passes: yes — 671/671 on Julia 1.10 LTS (unchanged from
  Chunk 17; this chunk touched only `docs/`).
- Docs build: `julia --project=docs docs/make.jl` succeeds locally,
  doctests pass.

## Next chunk

**Chunk 19 — Deferred-work inventory.** File GitHub issues for
everything out-of-scope: GUI, all variant scripts, low-RAM volume
variant, analysis-and-plotting helpers, experimental utilities, MEX
self-tests, plus the port-specific deferrals accumulated in chunk
notes. Each issue should reference the upstream `.m` files and the
algorithmic role.

## Watch out for

- **Chunk 19 produces GitHub issues, not code.** Check whether the
  repo has a GitHub remote and `gh` is authenticated before assuming
  issues can be filed; if not, the deliverable may need to be a
  `DEFERRED_WORK.md` inventory document instead. Confirm the intended
  output form.
- **The deferral list is already substantial** — scattered through
  every chunk's "Deviations from upstream" section in
  `ANALYSIS_PLAN.md` (e.g. temporal-focusing scattering, cortical-
  light-path orchestrator, `scan_ideal`/`single_scan_stack`,
  `make_avi`, `tifinitialize`/`tifappend` streaming, L1-penalised
  `times_from_profs`, polyphase resampling). Harvest these alongside
  the "Out of scope" section at the bottom of the plan.

## Working stance reminder

`ANALYSIS_PLAN.md` "Working stance" authorises autonomous chunk
progression with auto-commits.

## Suggested next workflow

`/clear` and re-run `/new-analysis-implement`.
