# Templates/LDMOS Stage-4 isothermal Id-Vd execution

Date: 2026-09-03  
Sentaurus: T-2022.03-SP2  
Scope: Phase A Stage 4, WP3-D decision, D5 classical isothermal entry

## Outcome

The complete Sentaurus D1--D5 single-factor chain ran successfully on the
sealed final SProcess mesh.  The hRecVelocity decision is closed as **do not
implement**: both gate curves are identical at the numerical noise floor and
there is no resolved source/drain component-current topology change.

The Vela D5 31-point qualification is not complete.  A contract-equivalent
gate prebias at `Vd=0 V` reached and saved both `Vg=4 V` and `Vg=8 V`, but the
first nonzero drain step at `Vg=4 V`, `Vd=0.01 V` failed the frozen block
absolute convergence gate.  Stage 4 therefore remains blocked and no Id-Vd
accuracy claim is made.

## Sentaurus single-factor results

All variants used the same `n1_fps.tdr`, `sdevice.par`, Fermi, OldSlotboom,
SRH/Auger, Math controls, contacts, two-dimensional width, voltage targets and
CurrentPlot grid.

| Parent -> child | Isolated factor | Vg=4 result | Vg=8 result | Decision |
| --- | --- | ---: | ---: | --- |
| D0 -> D1 | self-heating | endpoint 19.43% | endpoint 48.55% | Phase B thermal work required |
| D1 -> D2 | hRecVelocity | max 3.84e-12% | max 4.35e-12% | Do not implement |
| D2 -> D3 | hQP declaration | max 3.66e-12% | max 3.64e-12% | Inactive in shipped Solve blocks |
| D3 -> D4 | eQP declaration | max 4.41e-12% | max 3.20e-12% | Inactive in shipped Solve blocks |
| D4 -> D5 | IALMob | median/P95 59.84%/60.74% | median/P95 21.96%/46.47% | Implementation required after D5 closure |

The D1--D2 contact-component audit found a maximum absolute change of
`1.50162e-17 A/um`, only `4.48043e-14` of the same-point terminal current, and
zero resolved sign changes.  This closes both branches of the pre-registered
hRecVelocity decision rule.

The hQP/eQP rows require a narrow interpretation.  The official IdVd deck
declares both quantum models in `Physics`, but its `Solve` blocks contain only
Poisson, Electron, Hole and, in D0, Temperature.  They do not solve either
quantum-potential equation.  The observed zero delta therefore proves that
these declarations are inactive in this template path; it does not qualify a
future active DG implementation.

## Vela D5 execution and failure classification

The qualified G3 contract was retained: exact mesh, external LDMOS AverageBox
carrier couples, material-local Poisson charge volumes, Fermi/OldSlotboom,
SRH/Auger, HighFieldSaturation, no IALMob/QP/avalanche, no predictor and no
physical parameter adjustment.

Three entry paths were tested:

1. Reclosing an old `Vd=0.1 V` G3 state directly at zero drain failed and was
   rejected as a non-equivalent path.
2. Adding the WP1.5 drain branch guard rejected a same-bias reclose because of
   a `14.084 mV` contact-majority QF re-anchoring.  It was not used to mask the
   mismatch.
3. Reproducing the Sentaurus path from all-zero equilibrium, ramping the gate
   at `Vd=0 V` through 17 points from 0 to 8 V, passed and produced hashed
   Vg=4/8 checkpoints.  This is the accepted Stage-4 entry.

From the accepted Vg=4 checkpoint, the zero-drain point passed without Newton
updates.  A bounded no-guard `0 -> 0.01 V` experiment then failed after 40
iterations:

| Metric | Value |
| --- | ---: |
| initial scaled residual | 1.41421356 |
| best residual / iteration | 0.22519356 / 3 |
| final residual | 1.40853869 |
| final Poisson / electron / hole blocks | 18.129 / 20034.716 / 4.9e-5 |
| top Poisson node | 695 at (-9.976355, 10.119157) um |
| failure class | `block_absolute_convergence` |

The residual falls initially, then oscillates and grows in the electron block;
smaller startup drain steps do not close it.  This is a WP1.5 extension for a
high-gate, near-zero-drain operating point, not an hRecVelocity or calibration
problem.

## Status and next action

- Sentaurus Stage-4 oracle ablation: **pass**.
- WP3-D hRecVelocity decision: **closed, do not implement**.
- D5 Vela gate prebias and zero-drain checkpoint: **pass**.
- D5 Vela nonzero Id-Vd continuation: **blocked at 10 mV**.
- 31-point Stage-4 curve gate: **not run to completion**.
- IALMob: **mandatory follow-on**, but must not be enabled before D5 solver
  closure.

The next minimal task is an assembly/Jacobian audit of the Vg=4 V, Vd=10 mV
transition, comparing iteration 0, best iteration 3 and final iteration 40.
It should first inspect the electron-continuity row scaling and mobility-field
Jacobian around the failing source-side/high-doping basin and node 695.  Only
after this transition passes the frozen residual gates should the two Id-Vd
curves be resumed.

Machine-readable neutral results are stored in
`reference_tcad/templates_ldmos_sentaurus2022/stage4_decision_summary.json`.
Raw decks, logs, TDR/PLT, state CSVs and normalized tables remain under the
ignored `reference_staging/templates_ldmos_stage4_*_20260903/` directories.
