# LDMOS Stage 4: production precision-preserving QF path (2026-09-09)

## Result and scope

The production residual, transport-cell-vector HFS Jacobian feedback, and
contact-current postprocessing now recover cell gradients from QF reference
differences plus increment differences. They no longer first collapse each
reference and increment into one absolute double. Newton restart and diagnostic
packing also subtract the old and new references before adding the increment.

The previously stalled Vg=8 V checkpoint closes with one Newton update under the
original gates. The qualified parent, compensated translated parent, and its
saved/reloaded state close without an update. Full-curve qualification is a
separate requirement; the single-point results do not establish a 0–40 V pass.

Worktree: `.worktrees/templates-ldmos-phase-a`, branch
`codex/templates-ldmos-phase-a`, starting HEAD
`959c1c91513c62caeb700c2f71d8cef194879419`. Existing sweep step-control changes
and earlier HFS Jacobian changes are retained, and included in run provenance.

## Implementation

- `include/vela/equation/AssemblerUtils.h`: common referenced cell gradient and
  existing area-weighted edge recovery. Cell geometry and edge weighting remain
  the same.
- `src/equation/CoupledDDAssembler.cpp`: referenced fields in residual, local
  diagnostics, SG edge diagnostics, and Jacobian assembly; the mobility feedback
  chain rule uses the same cell gradient. Explicit frozen QF substitutions retain
  their specified physical-vector semantics.
- `src/post/ContactCurrent.cpp`: the same recovery from saved QF references and
  increments, with physical QF values used when no referenced representation is
  present.
- `src/solver/NewtonSolver.cpp`: restart and diagnostic repacking evaluate
  `(old_reference - new_reference) + increment` in extended precision. Adding
  the tiny increment to the large old reference first would lose it even before
  the subtraction in some checkpoints.
- `scripts/translate_dd_state.py`: reusable CSV frame translation with explicit
  transport-node ownership and compensated reference relocation. It preserves
  density/quantum columns, validates finite values and paired reference fields,
  and refuses output overwrite.

These changes address the classical Fermi–Dirac, vector-HFS LDMOS path used
below. They do not constitute a qualification of unrelated quantum or avalanche
models. No model defaults or acceptance tolerances are relaxed.

## Configuration and units

UCRT64 GCC, C++20, CMake `windows-ucrt64-release`; Eigen SparseLU/COLAMD is the
selected solver backend. Configuration detects HDF5/TDR and SuiteSparse
UMFPACK/SPQR support; their availability does not change the selected backend.

Mesh: 10,241 nodes, 19,782 triangles, 30,022 edges. Input geometry is in µm,
density in m^-3, voltage in V; internal unit scaling uses cm, and reported
terminal current is A/µm. Physics: 300 K, Fermi–Dirac statistics, Old Slotboom
BGN, SRH/Auger, `constant_field` mobility with `transport_cell_vector` QF-gradient
HFS, AverageBox/barycentric volumes, material-local Poisson, legacy node-local
contacts. Avalanche, quantum, and predictor are off.

Original gates are retained: Poisson/electron/hole block residual ceilings
5e-8 / 1e-11 / 3e-10; local carrier-row epsilon 1e-8 with zero violations;
terminal KCL ratio 1e-8. Curve children retain the 160-update budget and original
reference-point, frame-equivalence, configuration, and lineage audits.

## Numerical evidence

Evidence root: `reference_staging/templates_ldmos_frame_fix_20260909/`.

| Production solve | Newton updates | Poisson residual | Electron residual | Hole residual |
| --- | ---: | ---: | ---: | ---: |
| Qualified 26.6667 V parent, original frame | 0 | 4.976964e-9 | 9.224689e-12 | 1.540313e-24 |
| Compensated translated parent | 0 | 4.925889e-8 | 8.801839e-12 | 1.620242e-23 |
| Translated parent saved and reloaded | 0 | 4.925889e-8 | 8.801839e-12 | 1.620196e-23 |
| Previously stalled translated state | 1 | 4.103401e-8 | 3.830474e-12 | 1.120469e-23 |

All four solves report zero local carrier-row violations. Exact results,
configuration/seed/binary hashes and per-contact currents are in
`production_single_points.json`; these are same-bias diagnostic solves, not new
curve points. The generic `current_total_A_per_um` field can refer to a different
terminal after frame translation; compare named drain/source currents instead.

An exact rational oracle verifies all 11,446 active electron/hole QF values from
the translated real checkpoint. Density and quantum columns are unchanged.
181 translated increments differ by one last bit from the prior Decimal-80
diagnostic fixture: the exact oracle confirms all 181 new values are correctly
rounded binary half-way cases. The historical fixture remains unchanged.
See `translation_fixture_check.json` and `translation_tie_audit.json`.

## Regression qualification

The new vector-HFS test deliberately places the limiter transition at a
sub-ULP QF drop, and checks electron/hole residuals, Jacobian blocks,
postprocessed currents, and residual-derived currents without a unit floor.
It uses Fermi–Dirac statistics, matching the target transport path. A second
test checks exact retention of increments smaller than the extended-precision
ULP through a zero-update restart in both reference frames.

Initial test development exposed two invalid test assumptions: using the
default Boltzmann path for the Fermi–Dirac tiny-drop flux comparison, and
demanding cross-frame residual identity after independently rounded ohmic
boundary construction. The first was corrected by specifying the target
statistics; restart correctness is tested directly by exact retained increments.
Residual/Jacobian/postprocessing frame invariance remains covered by the first
test. Production code and acceptance tolerances were not changed for these
test corrections.

Seven Python translation tests pass. First full CTest run: 737/739 passed,
with only the two then-unadjusted new cases failing; initial logs are preserved
as `ctest_initial.log` and `ctest_initial.xml`. The corrected isolated C++ cases
pass. Final registered Release CTest: **739/739 passed**, 263.66 seconds
(`ctest.log`, `ctest.xml`). Both final new cases fail against the preserved
pre-fix core (27 assertions, 8 fail), confirming regression sensitivity;
see `negative_control.json` and `negative_tests.log`. `git diff --check` passes.

Frozen production runner SHA256:
`968ef00ed4cc16447c98853b4044f0895238f77f4972e18399694c7546cad97d`.
`environment.json`, `binaries.json`, and `source.diff` capture the exact source,
runner and successful test provenance before fresh continuation.

## Fresh curve qualification status

The real-mesh frame change at 26.6667 V has passed the original qualification:
normalized terminal-current difference 3.637536e-14; maximum electron/hole
density relative differences 2.132235e-13 / 2.132961e-13; KCL ratios
6.935784e-13 / 7.660597e-13. All potential/density/current/KCL equivalence
checks pass (`frame_to_28/fixed/ledger.json`, qualified frame event).

Continuation from 26.6667 to **28 V completed at 16:23:21** with 251 Newton
updates, zero rollbacks, 46 children, 23 transfers and both exact endpoints.
All 25 distinct accepted states pass the original qualification and immutable
source/binary/configuration/seed lineage audit (`frame_to_28/audit_summary.json`,
`integrity_pass: true`). This establishes successful continuation beyond the
previous stopping point, not merely a same-bias closure.

The Vg=4 0–1.3333 V smoke check completed at **16:29:50** with 214 Newton
updates, zero rollbacks, 45 children and both exact endpoints. All 24 distinct
accepted states pass the integrity audit
(`vg4_zero_to_1p333/audit_summary.json`).

The fresh **Vg=8 0–40 V** curve started at **16:29:50** from the original
qualified 0 V seed. All 31 exact reference points are required. Its live
provenance/status is in `queue_progress.json`, `vg8_full/progress.json`, and
`vg8_full/fixed/ledger.json`. Full-curve acceptance remains pending; this run
must not be combined with the older-binary Vg=4 full curve to claim a new
two-gate D5 pass.

At the user's request, a fixed deadline is **2026-09-09 16:43:36 Asia/Shanghai**,
30 minutes after the pause request was acknowledged. `run_with_deadline.py`
owns the validation/full-curve controller process tree and terminates that
tree at the deadline if unfinished, preserving accepted checkpoints and logs.
`deadline_progress.json` records the result. The existing `ldmos-d5` current-chat
heartbeat is active at a one-minute interval and verifies real process activity.
A deadline pause is not a numerical rejection and must not be reported as a
completed full curve. No automatic restart is authorized.
