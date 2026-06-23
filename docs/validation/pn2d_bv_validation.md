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
