# Templates/LDMOS G3 maximum-gm reclose and Jacobian audit (2026-08-31)

## Scope

This follow-up starts from the two Sentaurus endpoint states on the failed
`Vg=1.0--1.16666666666667 V`, `Vd=0.1 V` maximum-gm segment. It retains the
qualified G3 contract: `legacy_node_local` source-short reconstruction,
external AverageBox carrier couples, interior GradQF plus contact-cell
ElectricField HFS fallback, Fermi statistics, OldSlotboom, SRH and Auger.
No production default is changed.

Each Sentaurus state is first reclosed at the same bias with the full Vela
Newton system. A production first-step probe and a symmetric finite-difference
JVP check are then evaluated on the unreclosed Sentaurus state. The localized
JVP directions use the union of the endpoints of the 12 largest edge-flux
feedback differences at each bias: 25 nodes directly and an 83-node one-cell
ring.

## Same-bias reclose

| Vg (V) | Sentaurus Id (A/um) | Frozen Vela operator Id (A/um) | Reclosed Vela Id (A/um) | Newton iterations | Final scalar residual |
|---:|---:|---:|---:|---:|---:|
| 1.000000 | 1.254155163e-6 | 1.254408686e-6 | 1.021797549e-6 | 5 | 4.07930e-12 |
| 1.166667 | 2.827916110e-6 | 2.828487876e-6 | 2.213015335e-6 | 5 | 9.13155e-12 |

Both recloses reproduce the independently swept Vela endpoint currents. The
endpoint slope separation is therefore stable under a complete nonlinear
reclose, not an artifact of the frozen-state current extractor:

| Segment slope | A/(um V) | Ratio to Sentaurus |
|---|---:|---:|
| Sentaurus | 9.442565684e-6 | 1.000000000 |
| Vela operator on Sentaurus states | 9.444475139e-6 | 1.000202218 |
| Reclosed/self-consistent Vela | 7.147306711e-6 | 0.756924225 |

## First Newton response

The first production Newton step already predicts the direction and most of
the magnitude of the complete electron-QF reclose:

| Vg (V) | Full reclose phin L2 (V) | Coupled step/full L2 | Coupled cosine | Carrier-only step/full L2 | Carrier-only cosine |
|---:|---:|---:|---:|---:|---:|
| 1.000000 | 0.0812557 | 0.971208 | 0.994735 | 0.020982 | 0.032702 |
| 1.166667 | 0.2104514 | 0.939838 | 0.980335 | 0.008404 | -0.389048 |

The growing endpoint error is thus already present in the full coupled local
linearization; the remaining four iterations refine rather than create the
branch change. In contrast, the carrier-only solve at fixed Sentaurus psi is
nearly orthogonal to the required update. It reduces the immediate electron
residual, but supplies only `0.84--2.10%` of the final phin norm. The missing
SSS first-step response therefore uses the Poisson/electron-QF cross blocks.
It does not establish that Poisson feedback causes the branch difference,
because the complete Sentaurus state mixes a distinct psi/contact residual
support. The VSV counterfactual below isolates that question. The unreclosed
electron-continuity norm rises from `7.84346e-4` to `1.99452e-3` across the
segment, while the full phin correction norm rises by `2.590x`.

## Local residual versus nonlocal state response

The 25 direct edge-feedback nodes account for `88.7006%` and `78.9247%` of
the free-node electron-residual L2 at the two endpoints. They account for only
`11.6991%` and `11.8667%` of the complete phin-update L2. The two largest
continuity rows at both biases are nodes 4571 and 4541 near
`(-9.9792, 1.54688) um` and `(-9.92108, 1.54688) um`; their residual magnitude
grows from about `5.09e-4/4.74e-4` to `1.15e-3/1.07e-3`.

Conversely, the largest final phin shifts occur away from those residual
maxima. At `Vg=1.0 V` the maximum is `2.803 mV` near node 4314; at
`Vg=1.166667 V` it is `6.085 mV` near nodes 4085/4086. This is direct evidence
that the current difference is produced by a localized continuity mismatch
propagating through the global inverse carrier Jacobian, rather than by a
local HFS amplitude error at the final high-shift nodes.

## Jacobian verification

For the direct 25-node phin direction, the analytic-versus-symmetric-FD JVP
relative errors are `2.32e-11` and `3.29e-11`. For the 83-node one-ring patch,
they are `1.97e-10` and `3.93e-10`; four individual focus-node directions are
also below `3.13e-14`. The tested electron-QF Jacobian is therefore internally
consistent at both real endpoint states. This does not prove equivalence to
Sentaurus's private Jacobian, but it excludes an analytic-versus-primal Vela
derivative defect on the localized support.

## Fixed-psi carrier-block reclose

The follow-up audit removes Poisson feedback from the experiment. At each
endpoint, the `VSV` state retains Vela `psi` and `phip` but substitutes
Sentaurus `phin`. Only the carrier block is solved; `psi` remains fixed, every
accepted step uses damping 1, and no physical parameter changes.

| Vg (V) | Initial electron residual L2 | Residual after first step | Final residual L2 | Steps | Final Id / independent VVV |
|---:|---:|---:|---:|---:|---:|
| 1.000000 | 6.545894e-3 | 4.541772e-4 | 2.371216e-10 | 4 | 0.999999999963 |
| 1.166667 | 1.805441e-2 | 2.745615e-3 | 4.687803e-10 | 4 | 1.000000000078 |

The first carrier step has cosine similarity `0.999951` and `0.999706` with
the complete VSV-to-VVV `phin` change. It reduces free-silicon `phin` RMSE
from `1.08865` to `0.04212 mV` (`96.13%`) and from `2.81958` to
`0.24035 mV` (`91.48%`). After four steps, endpoint `phin` matches the
independent VVV state to RMSE below `7.8e-17 V`; `phip` changes by less than
`7.0e-17 V`. The resulting gm is `7.147306712e-6 A/(um V)`,
`1.00000000018x` the production Vela gm and `0.756924225x` Sentaurus gm.

At fixed Vela electrostatics, the Vela electron-continuity fixed point itself
therefore reproduces the full gm deficit. Poisson feedback, hole-QF feedback
and line-search branching are not required to create it.

## Spatial source of the growing phin response

The follow-up classifies nodes by exact mesh membership. At `VSV`, the
fraction of electron-residual L2 energy on the Si/dielectric interface is
`50.66%` and `60.23%`; including one silicon cell ring raises this to `91.92%`
and `95.20%`. Contact one-rings carry only `1.31%` and `0.87%`.

The largest rows at both endpoints form the gate-channel interface cluster
around nodes 3721, 3974, 3973, 4091, 3727, 4021 and 3771
(`x=-9.98 to -9.88 um`, `y=3.35--3.94 um`). Its `VSV` residual L2 grows
`2.75813x` across the failed gm segment. In contrast, the complete Sentaurus
state (`SSS`) retains a separate source/drain contact-neighborhood component:
contact-one-ring energy is `89.33%` and `70.24%`, and its first-step cosine
toward VVV is only `0.0327` and `-0.389`. It must not be used as the causal
support for the Vela branch shift.

The dominant `VSV` rows have diagonal fraction near `0.5` and
off-diagonal/diagonal magnitude near `1`, consistent with conservative
transport rows rather than an isolated singular row. The bias-growing source
is localized to the gate-channel interface-box continuity balance and its
global inverse response, not the source-short contact rows.

## Extended Jacobian verification

Three one-ring directions are tested on each `VSV` state: the dominant
interface `phin` patch, the separate SSS contact patch, and a coupled
`psi-minus-phin` interface direction. Analytic-versus-symmetric-FD relative
error is below `3.9e-9` for pure-`phin` directions and below `1.07e-6` for the
coupled direction. Together with the earlier edge-patch checks, this excludes
an implementation error in the tested mobility-field/carrier Jacobian.

## Decision and next target

1. The two-endpoint electron-QF reclose task passes. Full coupled reclose takes
   five iterations; fixed-psi carrier reclose takes four full steps and
   independently reproduces both swept Vela endpoints.
2. The maximum-gm deficit is classified as a Vela electron-continuity
   fixed-point difference concentrated in the Si/dielectric interface one-ring.
   It is not an AverageBox/HFS terminal operator-slope error, a Poisson or hole
   feedback effect, or an analytic-Jacobian inconsistency on the tested support.
3. The next minimal causal audit is coefficient-level comparison of the
   gate-channel interface-box electron rows: signed edge flux, control-volume
   weighting, source term and row scaling, using Sentaurus equation-balance
   evidence where available. The SSS contact-neighborhood residual remains a
   separately stratified effect.
4. The known-difference entry remains draft; no IALMob, predictor, mobility,
   threshold, HFS or contact parameter tuning is authorized.

## Reproducibility

- Driver: `scripts/audit_templates_ldmos_g3_gm_reclose.py`.
- Unit test: `tests/regression/test_audit_templates_ldmos_g3_gm_reclose.py`.
- Fixed-psi carrier driver:
  `scripts/audit_templates_ldmos_g3_phin_reclose.py`; unit test:
  `tests/regression/test_audit_templates_ldmos_g3_phin_reclose.py`.
- Runtime evidence:
  `reference_staging/templates_ldmos_g3_gm_reclose_jacobian_20260831/`.
- Same-bias reclose states:
  `reference_staging/templates_ldmos_g3_gm_reclose_20260831/`.
- Fixed-psi carrier/JVP evidence:
  `reference_staging/templates_ldmos_g3_phin_reclose_20260831/`.
- Frozen endpoint evidence:
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1_20260831/`
  and
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831/`.
