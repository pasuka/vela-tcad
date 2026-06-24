# PN2D BV Validation Methodology

## Source Of Truth

The source deck is
`reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd`. It uses
`Avalanche(VanOverstraeten)` and sweeps Anode to `-20.0 V`.

## Artifact Refresh

Sentaurus artifacts are refreshed through
`scripts/run_sentaurus_vm_reference.py pn2d --stages 0v,iv,bv`. Dry-run planning
uses `--dry-run`; live upload/run/fetch omits `--dry-run`.

## Required Artifact Checks

The BV refresh must contain `pn2d_bv.plt`, a clean BV log artifact such as
`pn2d_bv.log` or Sentaurus-generated `pn2d_bv.log_des.log`, and 201
endpoint-inclusive `pn2d_bv_multibias_*_des.tdr` files.

## Model Decision

When the refreshed source remains `Avalanche(VanOverstraeten)`, Vela's existing
`van_overstraeten` implementation is the target and Okuto-Crowell remains
contrast-only.

The refreshed Sentaurus BV source and currently available generated artifacts
use `Avalanche(VanOverstraeten)`. No Okuto-Crowell model is added in this
validation pass.

## Comparison Layers

1. Curve comparison checks current trend and documented windows.
2. Field comparison checks potential, electric field, carrier density, mobility,
   and avalanche generation/source density at selected biases.
3. Coefficient checks compare alpha to alpha; generation checks compare
   generation/source integral to generation/source integral.

## Accepted Gates

The promoted automated gates are VM-free and lightweight: parser/import
provenance, artifact validation with synthetic fixtures, BV max-field trend, and
documented comparison summaries.

## Non-Goals

This pass does not claim full `0..-20 V` absolute-current parity, does not
promote hidden scalar source calibration, does not rewrite SG flux divergence,
and does not add LDMOS/IGBT/MOS BV validation.

## High-Bias Interpretation

Windowed current agreement and high-bias knee shape are reported separately. If
knee-shape evidence remains divergent, the methodology records it as an open
physics/parity limit rather than hiding it behind a broad current band.

## Current Diagnostic Result

The derived Vela `-20.0 V` candidate deck under
`build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_validation_candidate`
currently records a controlled non-convergence diagnostic before the requested
high-bias comparison points. The latest local run reached a last stable bias of
about `-0.2301056002 V` and then failed at `-0.2301056003 V` with
`max_iterations`. The generated comparison report therefore treats the requested
`-0.5, -2, -5, -10, -20 V` candidate curve and field points as missing and
leaves the `-13.2..-13.0 V` current window `not_evaluated`; this is diagnostic
evidence, not a promoted `-20 V` parity gate.

After the live `pn2d_bv_validation_refresh` Sentaurus VM run, the refreshed BV
artifacts validate with 201 endpoint-inclusive multibias TDR files and a clean
BV log artifact. The local Vela candidate still does not reach `-20.0 V`.
Manual cross-checks against the Sentaurus Device user guide led to one code
alignment in this pass: GradQuasiFermi avalanche coefficient fields now fall
back to electric field in contact-touching elements. Control runs show that this
is not sufficient for the high-bias continuation limit: `electric_field` driving,
legacy `mobility_density_gradient`, `impact_ionization = none`, and a higher
Newton iteration budget all still stop around `-0.20..-0.23 V`. The remaining
gap is therefore treated as a coupled DD continuation/contact-state problem, not
as resolved BV physics parity.

A 2026-06-23 three-point localization (`c04edbf`, `51d8bf1`, current branch)
now narrows the remaining blocker to the `c04edbf -> 51d8bf1`
Newton/transport/globalization change set. Under `max_update=0` and
`quasi_fermi_update_limit_V=0.1`, the older baseline reaches `-3 V`, while the
current branch still fails before the canonical `-3 V` gate: impact-on fails at
`-0.40625000011175866 V` with `line_search_non_decrease`, residual
`1.9337e-8`, raw step norm `107.5935`, and positive finite carriers; the
no-impact control fails at `-0.8165000001609326 V` with `max_iterations`,
residual `8.5834e-9`, raw step norm `113.1072`, and positive finite carriers.
The failed Poisson-recorrection clamp trial was not retained. The follow-up
step probes show that the large step is carrier-block dominated:
Poisson-only steps are `~1e-8`, while carrier-only capped steps remain
`108..113` because hundreds of quasi-Fermi nodes hit the `0.1 V` cap after raw
carrier-row solves with near-zero diagonal dominance. Jacobian block finite
differences match the analytic blocks to about `1e-5` or better, so the next
work should target carrier-block conditioning/floor-row policy rather than
avalanche or Poisson recorrection.

The follow-up localization isolated the practical continuation blocker to the
high-field mobility field-sensitivity path in the Newton Jacobian. Keeping the
`masetti_field`/`quasi_fermi_gradient` residual physics but setting
`solver.mobility.jacobian_field_derivatives = false` unblocks the local
continuation gates: impact-on and no-impact both reach `-3 V`, and the impact-on
staged sweep reaches `-5 V`, `-10 V`, and `-20 V` with 401 accepted points and
3 Newton iterations at the final point. This updates the diagnostic status from
"cannot reach `-20 V`" to "numerically reaches `-20 V` with frozen high-field
mobility Jacobian." Full BV acceptance remains open after the acceptance refresh. The frozen-Jacobian
visual run converges `0, -0.5, -2, -5, -10, -20 V`, and the `-13.2..-13.0 V`
current-window gate passes with ratios around `0.804`. The high-bias gate still
fails: at `-20 V`, Vela current remains low by `0.8918` decades, the Vela curve
has no one-volt growth threshold above `1.5` or `2.0` in `-10..-20 V` while
Sentaurus reaches those thresholds at `-19 V` and `-20 V`, and field/state
parity still shows large high-bias density, electric-field, and thresholded
avalanche-generation errors. The next task is therefore high-bias avalanche and
carrier-feedback parity, not Newton reachability.

A follow-up frozen-Jacobian high-bias feedback localization keeps this status
unchanged but narrows the open physics axis. At `-20 V`, electric-field magnitude
is already close to Sentaurus, and the corrected continuity diagnostic uses the
actual material `ni = 1.4638914958767616e16 m^-3`; effective `ni` then matches
Sentaurus-inferred values at the active endpoints. The remaining current/source
deficit tracks absolute quasi-Fermi/carrier-density state offsets of roughly
`47..48 mV` for `psi-phin` and `56 mV` for `phip-psi`, causing about `0.8..0.95`
decades of carrier-density deficit. The next BV task is therefore a minimal
absolute-state feedback/branch-alignment probe, not alpha(E), material-ni, or
hidden source-scale calibration.

The minimal absolute-state feedback probe confirms this interpretation. On the
active `-20 V` edges, scaling Vela source density by the endpoint carrier-density
factors reconstructed from Sentaurus absolute `psi/phin/phip` state recovers
about `0.883` decades of source gap and moves the focus-edge source proxy from
`-0.840` to `+0.043` decades relative to Sentaurus. This keeps the next
investigation focused on high-bias absolute quasi-Fermi/carrier-density branch
selection.
The follow-up branch-offset probe further rules out a contact Dirichlet or global
electrostatic gauge explanation. Contact nodes remain aligned at `-20 V`, while
non-contact impact-active nodes have median `delta(psi-phin) = -0.04733 V` and
`delta(phip-psi) = -0.05547 V`, with electron/hole density deficits of about
`0.795/0.932` decades. The offset is negligible through `-5 V`, small at
`-10 V`, and large only by `-20 V`, so the next production-facing work should
probe the interior high-field carrier-continuity branch.

The active-support QF shift replay then tested that next hypothesis directly.
At the `-20 V` 99th-percentile avalanche support, Sentaurus and Vela active
nodes have zero overlap (`20` false-negative and `20` false-positive nodes), so
this replay is diagnostic rather than corrective. Globally shifting Vela's QF
branch by the measured active-state offsets moves Sentaurus-only active-node
transport ratios from `0.0955/0.1185` to `0.596/1.013` for electron/hole terms,
confirming branch-state causality. A hard support-only shift overshoots badly
(`14.20` hole transport on Sentaurus-only nodes and `8.424` electron transport
on Vela-only nodes), so the next production-facing task is a smooth active-region
or continuation-level branch-control experiment judged by the curve-level
`-10..-20 V` knee gate.
The smooth branch-control backscan tested the production-facing version of that
idea. A Gaussian QF shift around Sentaurus false-negative active support
(`decay_length_um = 0.05`) preserved contact nodes and converged a real DCSweep
backscan from `-20 V` to `-10 V`, but it did not move the knee: `-20 V` current
changed by only `0.00034` decades versus the frozen visual baseline and remained
`-0.8914` decades below Sentaurus. Residual probes show the shifted state is
Poisson dominated while carrier residuals at false-negative support remain tiny,
so the next task is a mixed-state Poisson/space-charge consistency audit of the
Sentaurus-like high-density active support.
The mixed-state charge audit rules out a direct Poisson/space-charge blocker for
that desired high-density support state. Replacing false-negative active-support
carrier densities with Sentaurus densities changes integrated net charge by only
`1.4116e-23 C/m` against `3.9116e-11 C/m` baseline net charge (`3.61e-13`). Even
replacing both false-negative and false-positive support changes only
`3.1196e-23 C/m`. The remaining branch issue should therefore be pursued in the
carrier-continuity flux/Jacobian balance, not by stronger QF shifts or Poisson
charge correction.

The follow-up active-edge replay and restart-relaxation probe narrows that
branch issue further. After fixing the SG-edge CSV loader to accept the current
`electron_flux_abs`/`hole_flux_abs`/`edge_area_m2` columns, the false-negative
active-edge replay at `-20 V` shows Vela baseline generation at only `0.138x`
Sentaurus, while a uniform Vela QF branch shift recovers it to `0.962x`; the
absolute QF density lever is therefore source-effective. However, a single-point
restart from the smooth shifted state converges in two Newton iterations back
toward the low-current branch: false-negative support retains only `0.198x` of
the electron QF shift and `0.255x` of the hole QF shift, with carrier-density
boosts reduced from `6.24x/8.55x` to `1.44x/1.73x`. The next experiment should
move coupled QF gradient/current-density branch state, not only the local
absolute density seed.
The coupled QF-gradient/current-density predictor experiment was then executed as
a stricter version of the branch-control probe. The new helper
`scripts/prepare_pn2d_bv_coupled_qf_predictor_state.py` starts from the Vela
`-20 V` state, selects active SG edges around requested support classes, blends
both endpoints' `phin/phip` to the Sentaurus endpoint QF pattern, and
reconstructs carriers from Vela inferred-ni while keeping `psi` fixed. The
mixed-state replay now accepts explicit `--state-csv-variant` inputs, so the
predictor initial state can be evaluated directly at the active-edge source
level.

For `false_negative` support only, the predictor touches `40` active edges and
`60` endpoint nodes. Its initial active-edge replay is source-effective:
false-negative generation is `1.0126x` Sentaurus and particle flux is `0.9826x`.
However the `-20 V` single-point restart again converges in two Newton iterations
to `-1.169034089095e-16 A/um`, effectively unchanged from the smooth-branch
restart. Endpoint relaxation shows the QF shift is mostly removed: false-negative
support retains only `0.2216x` electron shift and `0.2462x` hole shift, with
carrier boosts reduced from `7.136x/8.158x` to `1.546x/1.677x`.

For the all-support predictor (`false_negative + false_positive`), the helper
touches `66` active edges and `92` endpoint nodes. Its initial replay restores
both active classes near Sentaurus (`false_negative` generation `1.0126x`,
`false_positive` generation `1.0224x`), but the single-point restart still lands
at `-1.169034109498e-16 A/um`. The all-support endpoint relaxation is the same
pattern: retained QF shift remains only about `0.22..0.26x`. This rules out
"move both edge endpoints to the Sentaurus QF/current-density pattern" as a
sufficient branch-control predictor. The next useful probe is the first Newton
linear solve/update on this predictor state: identify which residual block or
Jacobian coupling drives the 75-80% QF-density rollback.

The first-Newton-step audit of the coupled predictor state was executed with
`scripts/diagnose_pn2d_bv_predictor_first_step_audit.py`. The helper converts the
predictor state to runner probe fields, runs `newton_step_probe` and
`newton_block_step_probe`, and reports rollback of the intended
`psi-phin`/`phip-psi` branch shift by support class.

For the false-negative-only predictor, the first full Newton step rolls back
`0.438x` of the electron branch shift and `0.418x` of the hole branch shift on
false-negative support. The carrier-only block step gives the same medians
(`0.438x` electron, `0.418x` hole), while the Poisson-only step is effectively
zero (`~2e-5` rollback fraction). The block residuals show why the solve still
accepts this direction: the predictor state is Poisson dominated
(`psi = 0.2443387`) but carrier residuals are already tiny (`phin = 1.41e-12`,
`phip = 2.01e-12`); the full step combines the Poisson correction with the
carrier block rollback.

The all-support predictor repeats the same pattern. False-negative support has
full-step rollback `0.438x/0.418x` for electron/hole, false-positive support has
`0.437x/0.425x`, and carrier-only matches full Newton. Poisson-only remains near
zero rollback. Thus the first-step rollback is not caused by Poisson/gauge or
Dirichlet movement; it is encoded in the carrier-continuity Newton block itself.
The remaining follow-up should inspect local carrier-row Jacobian coefficients
and RHS signs on the active endpoints, especially why a source-effective QF
branch is treated as a carrier-continuity residual reduction direction back
toward the low-density branch.

The active-endpoint carrier-row audit has now been executed with
`scripts/diagnose_pn2d_bv_predictor_carrier_row_audit.py`. The helper compares
the same support nodes in the baseline `-20 V` state, the coupled predictor
initial state, and the first Newton trial state by running
`newton_carrier_row_probe` and `newton_carrier_term_probe`, then joining row
stiffness, raw QF updates, flux, recombination, impact, and term-sum residuals.

For the false-negative-only predictor, the predictor state raises the median
false-negative electron flux term by `7.09e-14` and the hole flux term by
`9.68e-14`, while the impact terms compensate only `-8.24e-15` for each carrier.
The resulting carrier residual increases are therefore still flux dominated
(`6.31e-14` electron, `8.85e-14` hole), and the raw row updates roll back
`0.438x/0.418x` of the intended electron/hole branch shift. The first trial
state reduces both flux and impact, but remains above baseline with
`0.341x/0.336x` raw rollback.

For the all-support predictor, the same mechanism appears on both active
classes. On false-negative support, predictor-minus-baseline median flux deltas
are `2.20e-14` electron and `2.07e-14` hole, impact compensation is
`-8.24e-15`, and residual deltas remain positive (`1.32e-14` electron,
`1.27e-14` hole). On false-positive support, electron/hole flux deltas are
`1.88e-14` and `1.54e-14`, impact compensation is `-8.45e-15`, and residuals
remain positive (`1.04e-14` electron, `6.99e-15` hole). Thus the rollback is
not caused by missing local source effectiveness; it is the carrier-continuity
flux/source balance itself treating the Sentaurus-aligned branch as an
over-fluxed state.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/carrier_row_audit_false_negative_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/carrier_row_audit_all_support_blend1/`

The next diagnostic should perturb the carrier-row flux/source balance directly,
preferably by an impact-source scale/sign sensitivity or an equivalent frozen
residual/Jacobian probe, before designing another predictor. The question is no
longer whether the desired QF branch creates source; it is whether Vela's
continuity residual/Jacobian gives that source enough feedback against the SG
flux increase.
The impact-source feedback sensitivity was executed on the same carrier-row audit
artifacts with `--impact-scale` in
`scripts/diagnose_pn2d_bv_predictor_carrier_row_audit.py`. This is an analytic
row/term sensitivity over the existing `newton_carrier_term_probe` decomposition:
for each active endpoint it reports the impact multiplier needed to close
`term_sum + (scale - 1) * impact` and the adjusted residual at requested scales.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_scale_sensitivity_false_negative_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_scale_sensitivity_all_support_blend1/`

For the false-negative-only predictor, the median scale needed to close the
predictor carrier rows is `7.97x` for electrons and `10.26x` for holes; after
the first trial it remains `8.47x` and `10.43x`. Thus that predictor is not only
source-effective but also severely under-coupled in the local impact feedback
needed to balance the SG flux increase.

For the all-support predictor, the required scale drops sharply. On predictor
false-negative support the required medians are `2.31x` electron and `2.37x`
hole; on false-positive support they are `2.07x` electron and `1.71x` hole. The
first trial remains in the same range (`~1.61..2.16x`). This says the all-support
state has a much better flux/source balance, but Vela still needs roughly a
factor-of-two stronger effective impact feedback at the active carrier rows to
neutralize the rollback.

Next, do not tune the predictor further. The useful production-facing audit is
now the impact feedback path itself: compare Vela's impact source sign, scaling,
volume/edge support, carrier-current weighting, and Jacobian coupling with the
reference TCAD conventions before considering any source multiplier or
continuation constraint.

## Impact Feedback Semantic Audit

The production-facing impact feedback audit has now been executed. The local
reference-code check found no sign reversal in the basic generation semantics:
Genius DDM forms `GII = alpha_n * |Jn| / q + alpha_p * |Jp| / q`, distributes it
with directional current weights, and injects the same generated pair source into
electron and hole continuity rows over the finite-volume truncated partial
volume. Charon computes `alpha_n * |Je| + alpha_p * |Jh|` and subtracts avalanche
from total recombination, i.e. treats it as generation. Vela likewise subtracts
the SG edge-current avalanche source from both carrier residual rows. The local
`devsim` tree did not expose an equivalent built-in impact-ionization assembly to
compare.

The concrete all-support predictor comparison is now:

- Active-edge replay: predictor false-negative generation is `1.01257x`
  Sentaurus and false-positive generation is `1.02244x` Sentaurus on the selected
  active x-edges.
- Source geometry replay, after updating
  `scripts/diagnose_pn2d_bv_source_geometry.py` to accept current C++
  `sg_avalanche_edges.csv` columns (`edge_area_m2`, `electron_flux_abs`,
  `hole_flux_abs`, `source_integral`), reports active endpoint area fraction
  `0.5` for both support classes.
- Multiplying those two facts gives effective active-support feedback of
  `0.506x` on false-negative and `0.511x` on false-positive nodes, matching the
  carrier-row sensitivity that needed roughly `2.31x/2.37x` and `2.07x/1.71x`
  electron/hole impact feedback to close the rows.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/source_geometry_all_support_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/active_edge_replay_all_support_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_scale_sensitivity_all_support_blend1/`

Conclusion: the all-support predictor's local active-edge avalanche physics is
already Sentaurus-sized; the remaining factor-of-two is explained by effective
finite-volume/support feedback into the carrier rows, not by a raw avalanche
coefficient sign or local QF branch-strength error. Do not introduce a hidden
source multiplier. The next task should inspect whether Vela's endpoint
`0.5 * edge_area` carrier-row injection should be compared to Genius/Sentaurus
truncated partial-volume ownership on active edges, and separately whether the
SG edge-current avalanche Jacobian should include source derivatives rather than
omitting them in the current Newton path.
## Impact Feedback Ownership Policy Probe

A follow-up ownership summary has been executed with the new helper
`scripts/summarize_pn2d_bv_impact_feedback_ownership.py`. The helper joins three
already generated all-support predictor diagnostics:

- active-edge replay generation ratio versus Sentaurus,
- source-geometry active endpoint area fraction,
- carrier-row impact scale needed to close the local residual.

Artifact:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_feedback_ownership_all_support_blend1/`

Real `-20 V` all-support results:

| support class | active-edge generation / Sentaurus | active endpoint area fraction | endpoint feedback / Sentaurus | full active-edge feedback / Sentaurus | required e/h impact scale |
|---|---:|---:|---:|---:|---:|
| false_negative | `1.01257x` | `0.5` | `0.50628x` | `1.01257x` | `2.30979x / 2.36797x` |
| false_positive | `1.02244x` | `0.5` | `0.51122x` | `1.02244x` | `2.07195x / 1.71073x` |

This makes the factor-of-two failure mode concrete: the local active-edge source
strength is already approximately Sentaurus-sized, while endpoint-half ownership
leaves the carrier row with only about half of that active feedback. The product
`endpoint_feedback * required_scale` lands near unity (`1.17/1.20` for
false-negative and `1.06/0.87` for false-positive), so the previous `~2x`
impact-scale sensitivity is consistent with source-support ownership rather than
a raw ionization coefficient error.

Current C++ status: the SG edge-current avalanche residual path injects a source,
but the analytic Jacobian path still explicitly omits those nonlocal edge-source
derivatives. That is a separate Newton-coupling issue. The next production probe
should therefore be ordered as:

1. First test a focused source-ownership variant for SG edge-current avalanche,
   comparing Vela endpoint-half ownership against a directional/truncated-volume
   policy analogous to Genius/Sentaurus active-edge ownership. This must remain a
   gated experiment, not a hidden multiplier.
2. Then test source-derivative Jacobian completion for the SG edge-current path
   against finite-difference block probes.
3. Accept neither change without the curve-level `-10..-20 V` knee-shape gate and
   local carrier-row residual gate moving in the same direction.
## SG Edge-Box Source Volume Probe

Executed the focused production probe proposed above by adding an explicit gated
configuration knob, `impact_ionization.source_volume_policy`. The default remains
`edge_half_box`, preserving the previous SG edge-current avalanche source
support. The probe value `edge_box` changes only the SG edge source support from
`0.5 * h * edge.couple` to `1.0 * h * edge.couple`; it is not a hidden avalanche
coefficient multiplier and does not add the still-missing SG source-derivative
Jacobian terms.

Real `-20 V` all-support single-point deck:

- Baseline config: `single_m20_all_support_blend1/simulation_bv_coupled_qf_predictor_all_support_blend1_single_m20.json`
- Probe config: `single_m20_all_support_edge_box/simulation_bv_coupled_qf_predictor_all_support_edge_box_single_m20.json`
- Sentaurus reference current at `-20 V`: `-9.10455666344e-16 A`

| case | converged | Newton iterations | current_total_A_per_um | abs ratio vs baseline | abs ratio vs Sentaurus | decade error vs Sentaurus |
|---|---:|---:|---:|---:|---:|---:|
| endpoint-half baseline | `1` | `2` | `-1.16903410949798e-16` | `1.0000x` | `0.128401x` | `-0.891432` |
| `edge_box` probe | `1` | `2` | `-6.39910455999440e-16` | `5.47384x` | `0.702846x` | `-0.153140` |

The max field is effectively unchanged (`560748.547233267` to
`560748.547233097 V/cm`), so the terminal-current movement is carrier-continuity
feedback from the source ownership policy rather than a field-state movement. The
probe closes about `0.738` decades of the `0.891` decade gap and leaves a
remaining Sentaurus multiplier of `1.42279x`.

Conclusion: source ownership is now confirmed as a production-relevant direction,
but the `edge_box` probe alone is not yet an acceptance change. The next ordered
work is: run the same gated policy through the `-10..-20 V` knee-shape gate and
local carrier-row residual gate, then independently test SG edge-current
avalanche source-derivative Jacobian completion against finite-difference block
probes.

## SG Edge-Box Backscan Knee And Carrier-Row Gate

The `source_volume_policy=edge_box` experiment was extended from a single `-20 V`
point to a matched `-20 -> -10 V` all-support predictor backscan. Two decks were
generated from the same all-support restart family:

- endpoint-half baseline: `all_support_blend1_backscan/simulation_bv_coupled_qf_predictor_all_support_blend1_backscan_m20_to_m10.json`
- edge-box probe: `all_support_edge_box_backscan/simulation_bv_coupled_qf_predictor_all_support_edge_box_backscan_m20_to_m10.json`

Both backscans converged with `11` integer-bias points. The knee-shape gate now
reports:

| curve | first 1 V growth ratio > 1.5 | first 1 V growth ratio > 2.0 | max abs log10 current error |
|---|---:|---:|---:|
| Sentaurus | `-19.0 V` | `-20.0 V` | n/a |
| endpoint-half baseline backscan | none | none | `0.891432` decades |
| `edge_box` backscan | `-16.0 V` | `-19.0 V` | `0.519580` decades |

The probe therefore moves the curve-level knee in the right direction but does
not pass acceptance: the `>2.0` growth point is still one volt early, and the
`>1.5` threshold appears too early at `-16 V`. Pointwise, `edge_box` is close at
some mid-window biases but overshoots around `-19 V` (`3.308x` Sentaurus) while
remaining low at `-20 V` (`0.703x` Sentaurus).

A local carrier-row audit was also run under the `edge_box` policy with the
matched `-20 V` states:

- `carrier_row_audit_all_support_edge_box_policy/`

On the active support classes, the endpoint-half state evaluated with `edge_box`
already has required impact scales below unity (`false_negative` e/h
`0.8898/0.8774`, `false_positive` e/h `0.7784/0.7033`). The converged `edge_box`
state pushes the same rows further source-strong (`false_negative` e/h
`0.6452/0.5778`, `false_positive` e/h `0.6568/0.4142`) and lowers carrier raw
rollback magnitudes. This confirms source ownership is a real lever, but the
full-edge policy is too coarse as a production acceptance change.

Next ordered task: introduce and test an intermediate/truncated ownership factor
rather than binary `0.5` versus `1.0`, or move to the independent SG source-
derivative Jacobian probe if the goal is Newton coupling rather than curve-shape
calibration. Any intermediate policy must be gated by the same backscan knee
summary and carrier-row required-scale summary.
## Intermediate Source-Volume Factor Probe

Implemented a default-off diagnostic override, `impact_ionization.source_volume_factor`,
for the SG edge-current avalanche source support. The value `0` preserves the
named `source_volume_policy` preset. A finite value in `[0.5, 1.0]` overrides
only the edge source support factor used in `factor * h * edge.couple`; it does
not change the Van Overstraeten coefficients and does not add the missing
source-derivative Jacobian terms.

The all-support predictor family was rescanned from `-20 V` to `-10 V` with
intermediate factors. Every single-point restart and 11-point backscan below
converged.

| curve | source-volume factor | first 1 V growth ratio > 1.5 | first 1 V growth ratio > 2.0 | max abs log10 current error |
|---|---:|---:|---:|---:|
| Sentaurus | n/a | `-19.0 V` | `-20.0 V` | n/a |
| endpoint-half baseline | `0.5` | none | none | `0.891432` decades |
| factor probe | `0.75` | none | none | `0.649696` decades |
| factor probe | `0.875` | none | none | `0.453697` decades |
| factor probe | `0.90625` | `-15.0 V` | none | `0.381288` decades |
| factor probe | `0.9375` | `-15.0 V` | none | `0.319208` decades |
| `edge_box` policy | `1.0` | `-16.0 V` | `-19.0 V` | `0.519580` decades |

The factor scan confirms that source support magnitude is a real lever, but it
also shows a branch/continuation effect rather than a smooth calibration axis.
Factors up to `0.875` improve absolute current error but still have no 1 V knee
in the `-20..-10 V` gate. At `0.90625`, the `>1.5` growth marker jumps to
`-15 V` without producing the Sentaurus `>2.0` marker at `-20 V`. Therefore an
intermediate scalar source-volume factor should remain diagnostic-only and is
not an acceptance fix.

Carrier-row audit for `0.875` shows the local active-support required impact
scales are close to unity on the converged predictor state (`false_negative` e/h
`1.0055/0.9722`, `false_positive` e/h `1.0222/0.9366`), while the global curve
still lacks the Sentaurus knee. The next ordered task is the independent SG
edge-current avalanche source-derivative Jacobian probe against finite-difference
block checks; if that does not move the knee, the remaining discrepancy is likely
in branch selection/current-continuation support rather than local source volume
alone.
## SG Avalanche Source-Derivative Jacobian Probe

Executed the independent source-derivative Jacobian gate after the intermediate
source-volume factor scan. Current C++ already carries an SG edge-current
avalanche source derivative path in `CoupledDDAssembler::assembleJacobian`: the
edge source is finite-differenced with respect to endpoint `psi` and the local
quasi-Fermi stencil, then injected into the electron and hole continuity rows.
The probe therefore checks whether that implementation matches finite-difference
block behavior on both a synthetic high-field fixture and real `-20 V` PN2D BV
states.

Synthetic `pn2d_jacobian_block_audit` result:

| bias | block | analytic norm | FD norm | diff norm | rel diff |
|---:|---|---:|---:|---:|---:|
| `-20 V` | `sg_avalanche` | `7.99055e18` | `7.99054e18` | `1.65184e13` | `2.06725e-6` |
| `-19 V` | `sg_avalanche` | `7.99055e18` | `7.99054e18` | `1.65185e13` | `2.06726e-6` |
| `-15 V` | `sg_avalanche` | `7.99055e18` | `7.99054e18` | `1.65170e13` | `2.06707e-6` |

Real all-support `-20 V` state block probes used
`simulation_type=newton_jacobian_block_probe`, contact bias `Anode=-20 V`, and
`blocks=["sg_avalanche"]`:

| state | analytic norm | FD norm | diff norm | rel diff |
|---|---:|---:|---:|---:|
| endpoint-half baseline | `6.43224692543827e-14` | `6.43223981306022e-14` | `1.07307058520190e-18` | `1.07307058520190e-18` |
| factor `0.875` | `3.30373904700992e-13` | `3.30373765609820e-13` | `5.08483411648269e-18` | `5.08483411648269e-18` |

Conclusion: the SG avalanche source-derivative Jacobian block is not the current
BV knee blocker. It is numerically aligned with finite differences on the
synthetic gate and on the real high-field states tested here. The next ordered
work should target branch/continuation support: compare predictor direction,
state handoff, and current-growth branch selection around the sharp transition
between factor `0.875` (no knee) and factor `0.90625` (`>1.5` marker jumps to
`-15 V`).

## Charon-Alignment Minimal Probes

Executed three minimal Charon-alignment probes on the current PN2D Sentaurus2018
BV deck, with a baseline row for comparison. All runs used the rebuilt release
runner, `max_update=0`, `quasi_fermi_update_limit_V=0.1`, adaptive reverse-bias
steps, and a `-3 V` target. Artifacts were written under
`build-release/reference_tcad/pn2d_sentaurus2018/reports/charon_minimal_probes/`.

| variant | converged rows | deepest reverse bias | terminal failure |
|---|---:|---:|---|
| baseline | `5/6` | `-0.875 V` | `line_search_non_decrease` |
| `minimum_field_V_m=5.0e6` | `4/5` | `-0.75 V` | `max_iterations` |
| `driving_force=effective_field_parallel_j` | `5/6` | `-0.8125 V` | `max_iterations` |
| secant predictor + branch guard | `9/10` | `-0.326835876953125 V` | `electron_density_p95_jump_exceeded` |

The probes are negative as production fixes. The Charon-style minimum-field
cutoff and current-aligned driving force do not recover the missing continuation
reach; both stop earlier than the baseline in this low-bias robustness window.
The branch-guard probe does catch a real electron-density branch jump
(`p95_abs_dex=3.3721` on the rejected row, threshold `2.0`), but it rejects much
earlier than the baseline and therefore acts as a diagnostic branch detector, not
as a continuation strategy.

Conclusion: the next useful task should focus on why the accepted Newton state
moves into that density branch before `-1 V`: inspect the first rejected/failed
transition row around `-0.3268..-0.875 V`, comparing predicted initial state,
accepted solution, and carrier-continuity residual direction. The Charon knobs
alone are not sufficient to reach the `-3 V` gate.
## Low-Bias Branch-Window QF-Cap Probe

Executed the follow-up branch-window diagnostic around the first rejected
transition from the Charon-alignment minimal probes. The branch-guard run was
repeated with `write_state_every_point_prefix` and Newton history enabled, then
the final two accepted states were used to reconstruct the secant predictor for
the rejected target bias.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_window_probe/branch_window_summary.json`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_window_probe/predicted_newton_step.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_window_probe/hot_nodes.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_window_qflimit_probe/summary.csv`

The rejected transition is:

| previous accepted | current accepted | rejected target | secant ratio | rejection |
|---:|---:|---:|---:|---|
| `-0.326834453125 V` | `-0.326835876953125 V` | `-0.326837585546875 V` | `1.2` | `electron_density_p95_jump_exceeded` |

The reconstructed predictor already sits on the density-jump threshold:

| state comparison | p95 electron-density jump | max jump | max node |
|---|---:|---:|---:|
| predictor vs current accepted | `2.0159389305 dex` | `2.0159399605 dex` | `753` |
| first Newton trial vs current accepted | `3.3721311137 dex` | `3.6958660389 dex` | `517` |
| first Newton trial vs predictor | `1.6799260723 dex` | `1.6799260785 dex` | `581` |

The amplified nodes are all in the p-type region (`net_doping=-1e23 m^-3` in
the sampled hot nodes). Their first Newton step has nearly zero `delta_psi` but
`delta_phin=+0.1 V`, so `delta(psi-phin)=-0.1 V`. This exactly corresponds to
`0.1/Vt/log(10)=1.6799 dex`, matching the additional density drop from predictor
to trial. The local `phin_residual` on these nodes is already tiny (`~1e-16`),
so the branch jump is not caused by a large local carrier-continuity residual;
it is the configured quasi-Fermi update cap being applied coherently across the
p-region minority-electron branch while the global residual decreases.

A minimal QF-cap sweep then tested whether this mechanism is actionable under
the same secant predictor and branch guard:

| `quasi_fermi_update_limit_V` | converged rows | deepest reverse bias | failure | max accepted p95 jump |
|---:|---:|---:|---|---:|
| `0.1` | `9/10` | `-0.326835876953125 V` | `electron_density_p95_jump_exceeded` | n/a |
| `0.05` | `88/88` | `-3.0 V` | none | `0.9832502361 dex` |
| `0.025` | `54/55` | `-0.36078397757152 V` | `line_search_non_decrease` | `1.9923007654 dex` |

Conclusion: the immediate low-bias blocker is a narrow continuation/update-cap
interaction, not the Charon minimum-field or current-aligned driving-force knob.
`quasi_fermi_update_limit_V=0.05` is the next best candidate to extend beyond
the `-3 V` gate, while `0.1` oversteps the minority-electron branch and `0.025`
becomes too restrictive for line search in this window.
## QF-Cap 0.05 Full-Window BV Probe

Extended the actionable low-bias candidate, `quasi_fermi_update_limit_V=0.05`,
from the `-3 V` branch-window gate to the full `-20 V` PN2D Sentaurus2018 BV
window. The run kept the same secant predictor and branch guard used in the
low-bias probe.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_minus20/qflim0p05_minus20.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_minus20/qflim0p05_minus20_last_state.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_minus20/knee_shape_minus20_to_minus10.json`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_minus20/summary.json`

Run result:

| setting | result |
|---|---:|
| target stop | `-20 V` |
| converged rows | `491/491` |
| deepest reverse bias | `-20 V` |
| runner exit code | `0` |
| last current | `-1.1683243405189043e-16 A/um` |
| max accepted p95 electron-density jump | `0.9832502361 dex` |
| max accepted electron-density jump | `1.1071203737 dex` |
| min terminal-current consistency ratio | `0.4046220518` |

The continuation robustness problem is therefore substantially improved by the
`0.05 V` QF cap: the same branch guard that rejected the `0.1 V` cap near
`-0.3268 V` now allows a complete `-20 V` sweep.

However, the curve shape remains on the low-current/no-knee branch:

| curve | first 1V growth ratio > 1.5 | first 1V growth ratio > 2.0 | max abs log10 current error, `-20..-10 V` |
|---|---:|---:|---:|
| Sentaurus | `-19 V` | `-20 V` | n/a |
| Vela qf cap `0.05 V` | none | none | `0.891695` decades |

Integer-bias current ratios against Sentaurus worsen toward breakdown: the log10
current error is `-0.0924 dex` at `-10 V`, `-0.5485 dex` at `-19 V`, and
`-0.8917 dex` at `-20 V`. The final rows have current-jump ratios very close to
unity, so this run confirms that `qf_limit=0.05` is a continuation-stability fix,
not an avalanche-knee fix.

Next ordered task: keep `quasi_fermi_update_limit_V=0.05` as the stable
continuation baseline and rerun the source-ownership lever on top of it, starting
with `source_volume_factor=0.875` and `0.90625`. The acceptance question is
whether the previously source-effective factors can now express a Sentaurus-like
knee without reintroducing branch jumps.
## QF-Cap 0.05 Source-Volume Factor Scan

Executed the first source-ownership scan on top of the stable
`quasi_fermi_update_limit_V=0.05` continuation baseline, keeping the same secant
predictor and branch guard. Both source factors completed the full `-20 V`
window.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_source_factor_scan/summary.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_source_factor_scan/factor_0p875.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_source_factor_scan/factor_0p875_summary.json`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_source_factor_scan/factor_0p90625.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_source_factor_scan/factor_0p90625_summary.json`

| `source_volume_factor` | converged rows | deepest reverse bias | last current | max abs log10 current error, `-20..-10 V` | first 1V growth ratio > 1.5 | first 1V growth ratio > 2.0 | max accepted p95 jump |
|---:|---:|---:|---:|---:|---:|---:|---:|
| `0.875` | `491/491` | `-20 V` | `-3.2026773385317461e-16 A/um` | `0.453746` decades | none | none | `0.9832501955 dex` |
| `0.90625` | `491/491` | `-20 V` | `-3.78381968433441e-16 A/um` | `0.381328` decades | none | none | `0.9832504328 dex` |

Compared with the `qf_limit=0.05` default source factor (`0.891695` decades
max error), the source-volume factor is now a stable and effective magnitude
lever. It does not, however, recover the missing Sentaurus-like avalanche knee:
both candidates still have no `>1.5` or `>2.0` one-volt growth marker in the
`-20..-10 V` window.

Next ordered task: extend the same stabilized scan to the higher source
ownership candidates, especially `source_volume_factor=0.9375` and the existing
`edge_box` ownership mode, to determine whether the knee only appears once the
source magnitude is closer to the previous best-fit branch or whether ownership
geometry is the remaining missing lever.

## QF-Cap 0.05 High Source-Ownership Boundary Probe

Extended the stabilized `quasi_fermi_update_limit_V=0.05` scan to the higher
source-ownership candidates requested after the initial `0.875` and `0.90625`
runs. The probe kept the same secant predictor and carrier branch guard.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/summary.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/factor_0p921875.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/factor_0p921875_knee_shape.json`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/factor_0p9375.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/factor_0p9375_newton_failure_diagnostics.json`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/edge_box.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/edge_box_p95_2p1.csv`

| variant | result | deepest reverse bias | max abs log10 current error, `-20..-10 V` | first 1V growth ratio > 1.5 | first 1V growth ratio > 2.0 | failure |
|---|---:|---:|---:|---:|---:|---|
| `source_volume_factor=0.921875` | `491/491` | `-20 V` | `0.349196` decades | `-20 V` | none | none |
| `source_volume_factor=0.9375` | `12` stable rows | `-0.6897608 V` | n/a | n/a | n/a | `line_search_non_decrease` |
| `source_volume_policy=edge_box` | `24` stable rows | `-0.4148548 V` | n/a | n/a | n/a | `electron_density_p95_jump_exceeded` |
| `edge_box`, p95 guard relaxed to `2.1 dex` | `30` stable rows | `-0.4469347 V` | n/a | n/a | n/a | `line_search_non_decrease` |

This places the stable high-source boundary between `0.921875` and `0.9375` for
this continuation setup. The `0.921875` run is the best stable scalar-source
candidate so far: it improves the max current error from `0.381328` decades at
`0.90625` to `0.349196` decades and introduces a late `>1.5` one-volt growth
marker at `-20 V`. It still does not reproduce Sentaurus, whose reference window
has `>1.5` at `-19 V` and `>2.0` at `-20 V`.

The full `edge_box` geometry is not just a branch-threshold edge case. With the
original `2.0 dex` p95 guard it rejects at `-0.4148566 V`; relaxing the guard to
`2.1 dex` only reaches `-0.4469361 V` before a Newton line-search failure. Thus
`edge_box` is too aggressive for the current low-bias branch support under
`qf_limit=0.05`.

Next ordered task: binary-search the scalar factor interval `0.921875..0.9375`
with the same branch guard, then inspect the first failing Newton transition at
the unstable side. The acceptance question is whether a slightly higher stable
factor can recover the Sentaurus `>2.0` marker or whether the low-bias Newton
branch failure is the hard limit for source-ownership calibration.

## Source-Volume Factor Branch-Divergence Probe (Task 1)

Read-only diagnostic
`scripts/diagnose_pn2d_bv_factor_branch_divergence.py` (defaults reproduce the
run below):

```
python scripts/diagnose_pn2d_bv_factor_branch_divergence.py
```

It compares the converged `qf_limit=0.05` BV curves that differ only in
`impact_ionization.source_volume_factor`. The three candidate decks are
byte-for-byte identical except for that scalar (`quasi_fermi_update_limit_V=0.05`,
secant predictor, branch guard `p95<=2.0 dex`, `driving_force=quasi_fermi_gradient`,
`step=-0.25`), so the only varying axis is the source-volume factor.

| factor | reaches -20V | first >1.5 | first >2.0 | max abs log10 current error |
|---:|:---:|---:|---:|---:|
| Sentaurus | n/a | `-19.0 V` | `-20.0 V` | n/a |
| `0.875` | yes | none | none | `0.453746` |
| `0.90625` | yes | none | none | `0.381328` |
| `0.921875` | yes | `-20.0 V` | none | `0.349196` |

Per-pair result over the `-20..-10 V` window:

| pair | branch diagnostics identical | branch divergence onset | current-magnitude divergence onset | max abs current ratio | classification |
|---|:---:|:---:|:---:|---:|---|
| `0.875` vs `0.90625` | yes | none | `-19 V` | `0.072417` dex | source-volume magnitude |
| `0.90625` vs `0.921875` | yes | none | `-14 V` | `0.146448` dex | source-volume magnitude |

Across the window every candidate shares an identical continuation branch:
same `predicted_initial_state`, same `branch_acceptance_status`, identical
`electron_density_jump_p95_abs_dex` trajectory (max `0.98325 dex` for all three),
and identical adaptive `retry_count`. The factors differ only by a smooth current
magnitude scaling. The `>1.5` marker that first appears at `0.921875` is therefore
a smooth threshold crossing of the rescaled magnitude at the window edge
(`-20 V`), not a discrete branch jump.

**Classification: source-volume MAGNITUDE on a single shared branch, not a
continuation/predictor branch selection.** Scalar source-volume selection cannot
synthesize the Sentaurus knee on this axis; the knee blocker lies in the
continuation/branch mechanism itself. The ordered next work is pseudo-arclength
continuation plus local minority quasi-Fermi update caps, not further scalar
factor binary search.

### Reproducibility correction

This probe does not reproduce the earlier "Intermediate Source-Volume Factor
Probe" table that reported `0.90625` giving a `>1.5` marker at `-15 V` and
`0.9375` reaching `-20 V` with a `>1.5` marker. Against the current on-disk
artifacts the on-disk knee summaries instead show `0.90625` with no `>1.5`
marker (`first_growth_over_1p5_V=null`), `0.921875` as the first factor with a
`>1.5` marker (at `-20 V`), and `0.9375` failing at low bias (`-0.6897606 V`,
`line_search_non_decrease`) so it never enters the `-20..-10 V` window. The
stale table should be re-derived from the reproducible scan before being cited.

## Minority Quasi-Fermi Update-Cap Scan (Task 3)

Harness `scripts/run_pn2d_bv_minority_qf_cap_scan.py` clones the BV deck with
default half-box avalanche source support (`source_volume_factor` removed) and
varies only the Newton quasi-Fermi update caps:
`solver.quasi_fermi_update_limit_V` (global/majority) and the new
`solver.quasi_fermi_update_limit_minority_V` (tighter, minority-only,
classified per node by net doping). All other continuation settings match the
`qf_limit=0.05` baseline (secant predictor, branch guard `p95<=2.0 dex`,
`driving_force=quasi_fermi_gradient`, `step=-0.25`).

Low-bias gate (`stop=-1 V`):

| global | minority | deepest reverse bias | gate result | max accepted p95 jump |
|---:|---:|---:|---|---:|
| `0.1` | uniform | `-0.327 V` | fail `electron_density_p95_jump_exceeded` | `1.9181` |
| `0.1` | `0.05` | `-1.0 V` | pass | `0.9833` |
| `0.1` | `0.025` | `-1.0 V` | pass | `0.7063` |
| `0.05` | uniform | `-1.0 V` | pass | `0.9833` |

Knee window (full `stop=-20 V`):

| global | minority | converged rows | reaches -20V | max accepted p95 jump | max abs log10 current error |
|---:|---:|---:|:---:|---:|---:|
| `0.05` | uniform (baseline) | `491` | yes | `0.98325` | `0.891695` |
| `0.1` | `0.05` | `491` | yes | `0.98325` | `0.891695` |
| `0.1` | `0.025` | `546` | yes | `0.70634` | `0.891695` |

**Result: improvement at the low-bias gate with no knee-window regression.** A
minority-only cap unblocks the exact `-0.327 V` gate that uniform `global=0.1`
fails, by suppressing the minority (electron) density jump that trips the branch
guard while leaving the majority carrier free. At full `-20 V` every minority
candidate reaches the same deepest bias as the uniform `0.05` baseline at an
identical `0.891695`-decade max current error; `minority=0.025` even lowers the
max accepted p95 jump (`0.706` vs `0.983`) and takes finer adaptive steps
(`546` vs `491` rows). The minority cap is therefore a continuation-stability
lever, not a knee-synthesis lever: like the baseline, none of these recover the
Sentaurus `>1.5`/`>2.0` markers (all `null`), consistent with the Task 1 finding
that the knee requires a different continuation mechanism (pseudo-arclength).

## Minority Quasi-Fermi Update Cap (Task 3)

The global `quasi_fermi_update_limit_V` clips electron and hole quasi-Fermi
Newton updates uniformly. The prior continuation sweeps showed this scalar is a
blunt instrument: `0.025` is too restrictive (low-bias `line_search_non_decrease`),
`0.1` fails near `-0.327 V` with `electron_density_p95_jump_exceeded`, and only
`0.05` reaches `-20 V`. The branch jumps that trip the guard are minority-carrier
quasi-Fermi excursions, so a uniform cap over-constrains the majority carrier.

A new optional Newton key `quasi_fermi_update_limit_minority_V` adds a per-node
minority-carrier cap. Each node is classified by local net doping
(`netDoping < 0` p-type makes `phin` the minority update; `netDoping > 0` n-type
makes `phip` the minority update). When set `> 0` it tightens only the minority
carrier to `min(quasi_fermi_update_limit_V, quasi_fermi_update_limit_minority_V)`
while the majority carrier keeps the looser global cap. The default `0`
reproduces the existing uniform behavior exactly (verified: the prior
`quasi-Fermi update limit recomputes Poisson correction` test is unchanged, and
a new `minority quasi-Fermi update limit caps only the minority carrier per node`
test confirms electron-minority clipping to the minority cap with the hole
majority retained at the global cap on a p-type node).

This is the localized continuation lever recommended after the Task 1 finding
that scalar source-volume selection cannot synthesize the knee. The end-to-end
BV comparison matrix (low-bias gate at global `0.1`/`0.05` with minority caps
vs. the uniform baselines, plus the `-20..-10 V` knee window) is the next
release-sweep experiment; the mechanism and its default-off safety are landed
and unit-tested here ahead of that campaign.

## Pseudo-Arclength Continuation (Task 2)

Tasks 1 and 3 localized the knee blocker to the continuation mechanism: at the
avalanche knee `dI/dV -> infinity`, so any voltage-parameterized step is
unstable there regardless of local source modeling or quasi-Fermi caps.
Pseudo-arclength (Keller) continuation treats the bias as an unknown and adds an
arclength constraint, keeping the augmented system nonsingular through the fold.

Landed and unit-tested in this change:

- **Generic engine** (`include/vela/simulation/PseudoArclength.h`): a tangent
  predictor plus a bordered-Newton corrector with the arclength constraint
  `N = x_dot . (x - x0) + theta^2 * lambda_dot * (lambda - lambda0) - Delta s`.
  The corrector reuses the system Jacobian through a `solveJacobian` callback
  (block elimination: `a = J^-1(-F)`, `z = J^-1(F_lambda)`,
  `Delta lambda = (-N - x_dot.a) / (theta^2 lambda_dot - x_dot.z)`,
  `Delta x = a - z Delta lambda`) with adaptive arclength stepping.
  `tests/test_pseudo_arclength.cpp` proves it traverses the unit-circle fold at
  `lambda = 1` onto the `x < 0` branch that voltage stepping cannot reach.
- **Device adapter** (`NewtonSolver::makeArclengthSystem`): builds an
  `ArclengthSystem` over the coupled drift-diffusion residual with the bias on a
  chosen contact as the continuation parameter, reusing the exact assembler,
  boundary-condition construction, and sparse Jacobian assembly of the Newton
  solve. `dF/dV` is a central finite difference of the contact boundary values.
  `tests/test_newton_solver.cpp` (`[newton][arclength]`) confirms the bordered
  corrector advances a converged PN-diode equilibrium off `0 V` while keeping
  `||F||_inf` below tolerance on the real device Jacobian.
- **Sweep config plumbing** (`sweep.continuation.arclength`, default disabled):
  parsed and bounds-validated in `DCSweep`, with `tests/test_dc_sweep.cpp`
  covering valid parsing and each rejection path. Default decks are unchanged.

**Honest status / remaining gap.** The numerical engine, the device-level
bordered corrector, and the default-off configuration surface are complete and
unit-tested. What remains before the BV acceptance markers (`>1.5` 1V-growth at
`-19 V`, `>2.0` at `-20 V`) can be reported is wiring the arclength driver into
the `bv_reverse` main loop (replacing the voltage-step iterator, feeding the
branch-guard p95 jump as a step-shrink signal, and emitting the existing CSV
columns) and running the `-20..-10 V` release sweep. That integration and its
multi-minute validation are deliberately staged after the corrector was proven
on the real device Jacobian here, to avoid destabilizing the shared Newton path;
the Sentaurus knee recovery itself is therefore not yet demonstrated and is the
explicit next step.
