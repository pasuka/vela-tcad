# SLOT observed element avalanche source integration — 2026-09-23

## Status and scope

Worktree: `slot-ldmos-bvds-next`; base commit `863f4ee7`. This report extends
the [contact kernel correction](slot_ldmos_contact_kernel_fix_2026-09-23.md).
The near-angle ambiguity is removed **for this fixed mesh** by independent
native observations. Remaining ppm low-field mobility differences are explained
by an effective screening-minimum clipping difference at 300 K. The standalone
kernel is connected to an explicit controlled source profile, including live
residual, local AD, source-only audits, global closure and process diagnostics.
Neither a universal native selector algorithm nor strict original BVDS is claimed.

All device experiments use Windows UCRT64 Release, UMFPACK,
`scaling.mode=unit_scaling`, original 13,881-node SLOT mesh (24,272 Si triangles),
300 K, restored doped PolySilicon, aligned material reference, Si-only SRH/Auger
with density enhancement, and the existing projected element IALMob transport.
Avalanche source geometry is separate from transport geometry. Original Newton,
all-active carrier-row and global closure thresholds remain unchanged.

## Geometry identification

Two artificial frozen states add distinct common potential/QF polynomials to
the original native 10 V state. They are solely PMI geometry probes: no Newton
solve, no physical checkpoint, and no continuation seed. Matching uses the
actual native output state after native contact projection.

Both states and both carriers independently identify the same vertex for all
24,272 Si elements. The original 10 V observations overlap on 19,519 elements,
with zero disagreement. Cell 9170 requires a saved-potential ULP error bound;
its field error is smaller than that bound while the alternative cell's error
is many orders larger. Cells 5618/5620 share a coordinate and near-identical
field; their distinct callback current vectors disambiguate them in both runs.
No per-element tie tolerance was fitted into the production kernel.

Evidence under `build/`:

- `slot_pmi_geometry_20260923/` and `slot_pmi_geometry2_20260923/`: native frozen
  runs, returned data and hash checks, both exit 0.
- `slot_geometry_choice_20260923/complete_selector_audit.json` and
  `finalize_geometry.py`: complete coverage and exceptional-match evidence.
- `slot_observed_assembly_20260923/native_avalanche_geometry.json`: production
  controlled-input copy, with mesh coordinates added to verify geometry identity.

## Contact mobility discrepancy

The native HFS-off frozen control preserves primary state and the original
PMI callback columns. A second PMI probe reads n, p, donors, acceptors, nearest
interface distance/normal/orientation. Densities and doping agree exactly;
all orientations are 100 and distance differences are only export roundoff.

Inverting a bulk-dominated corner mobility gives effective native G floors:

| Carrier | Vela mathematical minimum | Inferred native effective floor | Upper crossing P |
| --- | ---: | ---: | ---: |
| Electron | 0.09810869654913681 | 0.09811015773601343 | 0.3258970688876528 |
| Hole | 0.0799431099966486 | 0.07994397371813261 | 0.2900551159556467 |

A separate held-out hole corner reproduces the inferred floor within about
3e-16. Applying those floors only in the independent diagnostic replay reduces
maximum positive-node HFS-off source differences from 2.895/6.539 ppm to
3.70e-12/1.07e-12 (electron/hole). Weighted L1 errors become 1.00e-15/3.03e-15.
This explains the observed numerical discrepancy; it does not identify the
native minimization algorithm for arbitrary temperature/mass. These inferred
constants are **not** installed in production.

Evidence: `slot_pmi_inputs_20260923/input_comparison.json`,
`slot_pmi_lowfield_20260923/inferred_screening.json`,
`screening_floor_upper_replay.py` and `screening_floor_upper_replay.json`.

Actual Release kernel replay with the observed selector and unchanged accurate
screening minimum yields:

| Frozen native bias | Electron weighted L1 source error | Hole weighted L1 source error |
| --- | ---: | ---: |
| 10 V | 2.3875e-14 | 1.5195e-12 |
| 37 V | 6.3287e-15 | 7.1438e-14 |

These compare local generation on a **frozen native state**, not terminal
currents or coupled device qualification. All 24,272 elements are included;
all nine-column source derivatives are finite. See
`slot_geometry_choice_20260923/release_avalanche_kernel_replay.json`.

## Controlled assembly and validation

Select `ialmob_element_density_grad_qf` explicitly; defaults and the existing
SG/GSS-Laux profile are unchanged. See [configuration schema](../config_schema.md)
for the atomic profile and supported parameters. The implementation checks
mesh identity/coverage and source measures; no saved mobility or density replaces
the live coupled state. Si-only source support also excludes restored poly.
Configuration fingerprints include live model parameters and prepared source
geometry. Alpha parameter-set/scaling controls are applied before SI conversion.

Regression scope includes every assembled potential column on SI and TCAD
fixtures, ordinary/contact-interface geometry, source AD versus finite
differences, process-source conservation, constrained/non-Si row exclusion,
unsupported model rejection, and malformed geometry/fingerprint checks.
Release full build passed. The complete Release CTest suite passed **902/902**
in 196.03 s. A final targeted scaling-linearity check brings the observed-source
tests to **622 assertions in 2 passing cases**. Its zero derivative direction
uses an operand-based finite-difference subtraction-roundoff bound, instead of
normalizing roundoff by a zero column. Nonzero derivatives retain the 3e-5
relative check. Evidence: `slot_observed_assembly_20260923/release_build.log`,
`ctest_release.log` and `observed_tests.log`.

Actual 10 V assembled audits include 24 initial directions and 12 adaptive
directions. Fixed macro perturbations cross QF-gradient norm branches when the
local QF span is extremely small. Branch-appropriate differences reduce every
tested carrier-block relative error below 1e-5. Tiny Poisson charge cross-columns
remain underresolved by double-precision differences; the full-column small
fixtures cover these analytically. This is not an exhaustive full-device
all-column finite-difference qualification.

## Qualified 2 V controlled-source checkpoint

The local aligned no-avalanche 2 V seed takes one complete Newton update, then
one carrier-only correction. The resulting state passes a **fresh full coupled
evaluation** with `initial_abstol` (zero further Newton updates), followed by an
independent solved-equation row audit. This is not an external-circuit reused
state or a cold-origin continuous scan.

| Check | Result |
| --- | ---: |
| Id | 9.264787215720655e-16 A/um |
| Existing native self-consistent avalanche reference | 9.2685914173646e-16 A/um |
| Current difference | -0.0410440% |
| Poisson block norm | 3.5936762402763475e-10 |
| Active electron / hole rows | 12634 / 12634 |
| Violating electron / hole rows | 0 / 0 |
| Independent maximum electron / hole row ratio | 1.13282e-6 / 5.39343e-5 |
| Electron / hole global closure ratio | 2.21899e-10 / 2.21899e-10 |

All reconstructed term sums exactly match reported row residuals. The reference
current is the previously recorded original native 2 V avalanche control in
[energy-reference validation](slot_ldmos_energy_reference_and_source_volume_2026-09-21.md).
`slot_observed_assembly_20260923/qualified_2v.json` records the checks and SHA-256
of runner, configuration, mesh, materials, doping, both geometries and state.

Direct 2-to-3 V and 2-to-2.125 V attempts each reached their 300 s execution
budget before producing complete I-V output. These are **budget interruptions**,
not network interruptions, proved numerical failures, or accepted checkpoints.
See `continue_from2/result.json` and `continue_small_steps/result.json`.

A subsequent bounded **24-update** attempt writes a complete rejected state
before external correction; the iteration budget changes, not physical or
acceptance thresholds. At 2.125 V its Poisson norm is 3.49087e-10 and global
closure passes, but 548 carrier rows fail (maximum ratio 1.78183). Four physical
carrier Newton corrections leave 466 violating hole rows (maximum 1.76956),
while all electron rows pass. Thus this is a demonstrated numerical continuation
failure, distinct from the earlier timed-out attempts.

The independent `carrier_regularization_scale=0` control fails after 21 Newton
updates, with 6,165 row violations and electron/hole closure ratios 1.914/1.871.
Four subsequent corrections also fail. It is not adopted. A density-coordinate
candidate from the original regularized failure is rejected with
`carrier_density has no residual-decreasing candidate`; no state is qualified.
Evidence: `continue_bounded24/result.json`, `regularization_off_control/result.json`
and `density_correct_2p125/local_corrector/step_001/stderr.log`, all under
`slot_observed_assembly_20260923/`. The last qualified new-source point remains
**2 V**. The follow-up must improve bias prediction/initialization and weak-hole
continuity before attempting further voltage continuation; increasing time alone
or disabling regularization has not established a solution.

## 10 V numerical limitation

`slot_observed_assembly_20260923/solve10/` failed numerically; subsequent
`correct10_without_poly/`, `prewarmed10/` and `correct10_original_seed/` attempts
also failed full coupled qualification. They completed without timeouts. In the
last original-seed attempt, both carriers have zero violating effective rows
(12,634 rows per carrier), maximum row ratio 2.399e-12, and global closure ratios
1.85e-16/9.25e-16. Yet the Poisson block norm is 1.27641e-9 and full Newton
reports `line_search_non_decrease`. Rejected I-V zero fields are placeholders,
not a zero-current device solution.

The limitation has two parts: the weighted Poisson residual is about 1.073e-9,
above `abstol=1e-9`, and the numerical-floor acceptance path is separately blocked
by the source-contact majority-hole QF drop 5.2492e-10 V (edge 3588--3589), above
its unchanged default limit 5e-11 V. Poly contacts are already near equilibrium.
Thus the floor guard cannot be replaced with a claim that carrier/global closure
alone establishes a converged device.

`roundoff10.py` makes bounded adjacent-representable-potential seed trials and
re-evaluates the **production** Poisson residual. The best probe norm falls only
from 1.27641e-9 to 1.25091e-9; subsequent full Newton plus all-row revalidation
still fails. This is a diagnostic seed experiment, not a tolerance adjustment
or accepted checkpoint. Its evidence is in `roundoff10/` and `roundoff10.log`.

Remaining device work is to separate electrostatic representation/roundoff
and nonlinear update limitations with independently checked arithmetic and a consistent full Jacobian,
then rerun the unchanged device gates. A higher-bias or load-line checkpoint
from the earlier projected-SG source cannot qualify this new source profile.
Full original-model BVDS and post-breakdown comparison remain unfinished.

## Resume record

At the end of this run all owned local solves, bounded controls and Release
checks have finished; no background continuation is left running. Changes are
in the local worktree. No production screening constant or tolerance was fitted
to the native result, and no failure/timeout file replaces the qualified 2 V
state. Resume from `qualified_2v.json`, checking its file hashes and retaining
the failed 2.125/10 V evidence for comparison. This run did not start additional
remote Sentaurus solves or restart the old paused monitor.
