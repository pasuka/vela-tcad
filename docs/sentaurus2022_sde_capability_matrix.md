# Sentaurus 2022 SDE/SDevice capability matrix

This matrix is normative for the paired SDE/SDevice workflow implemented by
`scripts/sentaurus_import.py device`. SDE semantics follow *Sentaurus Structure
Editor User Guide*, release T-2022.03. A feature marked unsupported must raise a
source-located error; it must never be silently approximated.

`offline-static` is the license-free default for public users, packaged
workflows and CI. No public acceptance path may require a Sentaurus executable,
license server or proprietary installation. `sentaurus2022-evaluated` is an
explicit opt-in reserved for internal development, reviewed oracle generation
and compatibility testing. Frozen, provenance-sealed oracle artifacts may be
consumed by license-free tests after review.

## Public Gmsh preprocessing milestones

`scripts/sentaurus_gmsh_pipeline.py` is the public, license-free SDE-to-mesh
path. Its default `--capability-level M2` is a frozen semantic subset, not a
claim that every capability in the broader paired workflow below is already
available through Gmsh.

| Capability | Gmsh M1 | Gmsh M2 (default) | Later milestone |
|---|---:|---:|---:|
| separated Rectangle + pick contact | yes | yes | - |
| per-body ABA/BAB and conformal interface groups | fail-closed | yes | - |
| straight-edge `insert-vertex` + complete split-edge contact | fail-closed | yes | - |
| ordered constant NoReplace/LocalReplace/Replace | NoReplace only | yes | - |
| directed line-baseline Gaussian + Eval/NoEval | fail-closed | yes | - |
| nonzero constant or analytical restriction decay | fail-closed | fail-closed | M3 |
| refinement placement/function | fail-closed | fail-closed | M3 |
| Polygon body, AB/ABiA/ABiB/XX, 3-D | fail-closed | fail-closed | excluded/unplanned |

For T-2022.03 constant placements, NoReplace is represented by omitting the
optional replacement argument. Replace or LocalReplace requires the numeric
decay-length positional slot first; use `0 "LocalReplace"` or `0 "Replace"`
when no decay is intended. The offline parser rejects the tempting but invalid
literal `"NoReplace"` form rather than accepting a deck that official SDE
rejects.

## Frontends and provenance

| Capability | `offline-static` (default) | `sentaurus2022-evaluated` (internal only) |
|---|---:|---:|
| Run `sde -e -l wrapper.scm` | no | yes |
| Enable `(set! evaluate-log-file #t)` before evaluating input | no | yes |
| Parse variable-expanded `.log.eval` | n/a | yes |
| Record source/evaluated SHA-256, SDE version, command, stdout/stderr and exit | source hash only | yes |
| Full Elk Scheme | no | delegated to installed SDE |
| `define` plus numeric expressions | yes | via evaluated log |
| `set!`, `if`/`cond`, `do`, `for-each`, `lambda`, `load` | fail-closed with internal-frontend hint | via evaluated log |

Both paths emit `vela.sentaurus_sde_frontend_manifest.v1`. T-2022.03-SP2 does
not flatten a nested `load` into `.log.eval`, so the evaluated path uses an
isolated wrapper containing the toggle followed by the unchanged, hash-sealed
source text. The evaluated path can additionally invoke `sde:build-mesh` and
preserve the official `_bnd.tdr`
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
| Positive/Negative/Both (`Symm`) | supported | T-2022.03-SP2 oracle: Positive is the left normal of the directed baseline |
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

Sentaurus Quasistationary `InitialStep`, `MinStep`, and `MaxStep` are normalized
by the original `Goal - start` span before Vela receives voltage steps. The
contract may shorten the stop bias without changing that normalization.
Terminal-current acceptance uses a mixed gate,
`abs(Ia + Ic) <= abs_limit + rel_limit * max(abs(Ia), abs(Ic))`, so roundoff at
near-zero current does not masquerade as a large physical conservation error.

Excluded: complete Elk reimplementation, ACIS, Sentaurus Mesh reimplementation,
3-D geometry, process emulation, particle profiles, and arbitrary external TDR
profiles.
