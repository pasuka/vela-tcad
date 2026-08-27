# SimpleMOS M8 original-physics SDevice validation (2026-08-27)

## Outcome

M8 execution is complete, but its original-deck comparison gate does not pass.
All 16 Sentaurus reference curves and all 16 Vela curves are qualified at the
exact 51-point gate lattice. Eight cases pass and eight fail the predeclared
numeric limits. The split is exact: all four `NWell=1e17 cm^-3` devices pass,
while all four `NWell=2e17 cm^-3` devices fail only at deep-off gate biases.

This result is an SDevice comparison. The seven additional process runs were
used only to materialize immutable upstream TDR inputs. No SProcess feature,
intermediate field, model, or result is validated here.

## Frozen contract

- Sentaurus Device release: `T-2022.03-SP2`.
- Source deck SHA-256:
  `3fd501cad368aa43c55a780253309428872a4b3349aa03727f8a875d04a4b536`.
- Physics: `OldSlotboom`, `PhuMob + HighFieldSaturation + Enormal`, and
  `SRH(DopingDependence)`.
- Biases: `Vd={0.05,1.0} V`, with exactly 51 direct `Vg=0:0.05:2.5 V`
  points and no interpolation.
- Numeric limits: maximum absolute log10 ratio `0.3 dex`, maximum relative
  error `1.0`, and matching trend; current floor `1e-18 A/um`.
- State policy: equilibrium -> drain ramp -> gate sweep, with each stage
  consuming only the preceding accepted state. Auxiliary drain points have no
  acceptance role.

The nominal n17 two-curve gate passed before any non-nominal Vela matrix run.
Its worst errors were `0.0422 dex / 10.20%` at `Vd=0.05 V` and
`0.0322 dex / 7.70%` at `Vd=1.0 V`.

## Matrix result

| Device | NWell (cm^-3) | GOxTime (min) | LDD dose (cm^-2) | Vd=0.05 V | Vd=1.0 V | Device result |
| --- | ---: | ---: | ---: | --- | --- | --- |
| n17 | 1e17 | 10 | 1e14 | pass: 0.0422 dex, 10.20% | pass: 0.0322 dex, 7.70% | pass |
| n18 | 1e17 | 10 | 2e14 | pass: 0.0421 dex, 10.17% | pass: 0.0329 dex, 7.88% | pass |
| n19 | 1e17 | 15 | 1e14 | pass: 0.0399 dex, 9.61% | pass: 0.0275 dex, 6.53% | pass |
| n20 | 1e17 | 15 | 2e14 | pass: 0.0457 dex, 11.10% | pass: 0.0280 dex, 6.66% | pass |
| n21 | 2e17 | 10 | 1e14 | fail: 0.3497 dex, 123.73% | fail: 0.4835 dex, 204.41% | fail |
| n22 | 2e17 | 10 | 2e14 | fail: 0.3211 dex, 109.45% | fail: 0.4470 dex, 179.87% | fail |
| n23 | 2e17 | 15 | 1e14 | fail: 0.5371 dex, 244.45% | fail: 0.5778 dex, 278.25% | fail |
| n24 | 2e17 | 15 | 2e14 | fail: 0.4732 dex, 197.27% | fail: 0.5510 dex, 255.63% | fail |

Matrix-wide qualification and comparison facts:

- 16/16 Sentaurus curves qualified, 51 direct points each (816 total).
- 16/16 Vela curves completed with all gate points converged.
- 8/16 comparisons pass; 8/16 fail.
- All 16 trends match and are nondecreasing.
- Worst observed error is `0.5778 dex / 278.25%` (`n23`, `Vd=1.0 V`).
- The maximum absolute `Vg=2.5 V` endpoint error is only `0.0132 dex`.

## Failure localization

Every threshold violation is confined to `Vg=0` or `Vg=0.05 V`:

| Case | Violating direct points | Violating Vg (V) |
| --- | ---: | --- |
| n21, Vd=0.05 V | 1 | 0 |
| n21, Vd=1.0 V | 2 | 0, 0.05 |
| n22, Vd=0.05 V | 1 | 0 |
| n22, Vd=1.0 V | 1 | 0 |
| n23, Vd=0.05 V | 2 | 0, 0.05 |
| n23, Vd=1.0 V | 2 | 0, 0.05 |
| n24, Vd=0.05 V | 2 | 0, 0.05 |
| n24, Vd=1.0 V | 2 | 0, 0.05 |

The failing currents are approximately `5e-17` through `1.1e-15 A/um`.
They are above the frozen `1e-18 A/um` comparison floor, so they cannot be
downgraded to diagnostic-only points. No limit or floor was changed after the
matrix result was observed.

The clean split at `NWell=2e17 cm^-3`, across both oxide times and both LDD
doses, is evidence for a concentration-sensitive deep-off parity gap. The
small strong-inversion and endpoint errors make the implemented
PhuMob/Enormal high-current path a less likely primary cause. This is a
localization, not yet a root-cause proof.

## Required follow-up before M9

M9 regression/documentation freeze is not ready. The next development task
should compare n20 and n21 at equilibrium and the failing gate points, with
special attention to:

1. `OldSlotboom` effective intrinsic density and band-gap-narrowing state at
   `NWell=2e17 cm^-3`;
2. minority-carrier densities, SRH generation/recombination, and contact
   equilibrium in the channel/well;
3. sub-femtoampere contact-current extraction and continuity residual balance;
4. whether the difference is already present before mobility field and
   Enormal terms become active.

Acceptance limits, the current floor, and the exact gate lattice must remain
unchanged during that diagnosis.

## Evidence

- Contract: `reference_tcad/simplemos_sentaurus2022/simplemos_m8_original_physics_contract_v1.json`
- Evidence manifest: `reference_tcad/simplemos_sentaurus2022/simplemos_m8_original_physics_evidence.json`
- Neutral references: `reference_tcad/simplemos_sentaurus2022/original_physics/`
- Comparison report: `reference_tcad/simplemos_sentaurus2022/original_physics/comparisons/comparison_report.json`
- Reproducible runner: `scripts/run_simplemos_m8_original_matrix.py`
- Regression test: `tests/regression/test_simplemos_m8_original_matrix.py`

TDR, PLT, Sentaurus logs, generated neutral exports, and Vela restart states
remain ignored and are not committed.
