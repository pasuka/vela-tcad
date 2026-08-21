# Sentaurus 2022 SDE/SDevice capability matrix

This matrix is normative for the paired SDE/SDevice workflow implemented by
`scripts/sentaurus_import.py device`. SDE semantics follow *Sentaurus Structure
Editor User Guide*, release T-2022.03. A feature marked unsupported must raise a
source-located error; it must never be silently approximated.

## Frontends and provenance

| Capability | `sentaurus2022-evaluated` | `offline-static` |
|---|---:|---:|
| Run `sde -e -l wrapper.scm` | yes | no |
| Enable `(set! evaluate-log-file #t)` before loading input | yes | no |
| Parse variable-expanded `.log.eval` | yes | n/a |
| Record source/evaluated SHA-256, SDE version, command, stdout/stderr and exit | yes | source hash only |
| Full Elk Scheme | delegated to installed SDE | no |
| `define` plus numeric expressions | via evaluated log | yes |
| `set!`, `if`/`cond`, `do`, `for-each`, `lambda`, `load` | via evaluated log | fail-closed with evaluated-frontend hint |

Both paths emit `vela.sentaurus_sde_frontend_manifest.v1`. The evaluated path
can additionally invoke `sde:build-mesh` and preserve the official `_bnd.tdr`
and `_msh.cmd` artifacts. It does not parse TDR as though it were Scheme.

## Geometry and contacts

| Capability | Status | Notes |
|---|---:|---|
| 2-D axis-aligned rectangle | supported | conformal tensor grid |
| ABA (New Replaces Old) | supported | mode captured on each newly created body |
| BAB (Old Replaces New) | supported | mode captured on each newly created body |
| AB, ABiA, ABiB, XX | fail-closed | never downgraded to ABA |
| partial/full containment and multiple old bodies | supported | zero-area remnants are removed from output regions |
| polygonal body, ACIS, 3-D geometry | fail-closed | outside current geometry model |
| `define-contact-set` display arguments | supported | preserved in IR |
| current/explicit contact selection | supported | missing explicit set is created and activated |
| `set-contact` edge and edge lists | supported | each pick must resolve to exactly one complete edge |
| `insert-vertex` on a straight outer 2-D edge | supported | splits the edge before contact expansion |
| endpoint, multi-edge, internal-interface ambiguity | fail-closed | no adjacent-cell majority heuristic |
| body contact, `remove`, 3-D `imprint` | fail-closed | diagnostic names the unsupported operation |

## Doping and mesh control

| Capability | Status | Notes |
|---|---:|---|
| constant profile to window/region/material | supported | source order and source location preserved |
| NoReplace | supported | accumulates with every active earlier profile |
| LocalReplace | supported | clears only the same species |
| Replace | supported | clears all earlier species at the sample point |
| constant Erf/Gaussian decay length | supported | boundary distance evaluated for line, rectangle and polygon windows |
| donor and acceptor aggregation | supported | kept separate until final CSV columns |
| Gaussian PeakPos + PeakVal | supported | with ValueAtDepth/Depth, Length, or StdDev |
| Gaussian Dose | supported | requires Length or StdDev |
| lateral Gauss/Erf | supported | Factor/Length/StdDev parameters preserved and evaluated |
| Positive/Negative/Both (`Symm`) | supported | Positive is the right side of the directed baseline |
| Eval/NoEval | supported | NoEval excludes the baseline itself |
| analytic region/material/eval-window restriction | supported | explicit restriction only |
| Line/Rectangle/Polygon Ref/Eval window | supported | polygon inclusion is boundary-inclusive |
| refinement region/material/window placement | supported | overlapping placement uses the smallest active spacing |
| max/min spacing | supported | min is a refinement stopping target, not a hard lower bound |
| MaxTransDiff/MaxGradient | supported conservatively | full placement refines to its stopping target, a safe superset of adaptive cells |
| MaxLenInt, DoubleSide, UseRegionNames | supported conservatively | full placement uses the requested interface length, including both sides |
| MaxInterval, multibox, offsetting | fail-closed | listed here and in source diagnostic |

## Unified execution and exclusions

`analysis.kind=iv` maps to `sweep.mode=iv`; `analysis.kind=bv` maps to
`sweep.mode=bv_reverse`. All parsed initial/sweep stages are copied into the
ordered `_sentaurus_execution_sequence`. Until cross-deck checkpoint handoff is
implemented, more than one Quasistationary sweep fails closed instead of losing
an earlier stage.

Run status is one of `generated`, `solved`, `accepted`, or `failed`. `accepted`
is written only after the checked-in `configs/sentaurus_e2e_contract.json` is
read and its mesh, curve, coverage, finiteness, monotonicity, and continuity
thresholds pass. Failure writes `fail_report.json`, removes stale success
reports, and records the last accepted bias.

Excluded: complete Elk reimplementation, ACIS, Sentaurus Mesh reimplementation,
3-D geometry, process emulation, particle profiles, and arbitrary external TDR
profiles.
