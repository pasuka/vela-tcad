# SimpleMOS SDevice validation plan

Date: 2026-08-26
Branch: `codex/simplemos-sdevice-validation`
Sentaurus release: `T-2022.03-SP2`

## Objective

Build a traceable Sentaurus Device versus Vela comparison for the official
`GettingStarted/swb/SimpleMOS` SDevice deck. Sentaurus Process output is an
immutable upstream input. This work does not implement or validate oxidation,
implantation, diffusion, deposition, etching, annealing, or any SVisual
extraction.

The target SDevice sequence is:

1. Poisson initialization;
2. coupled Poisson/electron/hole equilibrium;
3. quasistationary drain ramp to `Vd`;
4. quasistationary gate sweep from 0 V to 2.5 V;
5. drain total-current comparison in A/um.

## Frozen source evidence

| Item | Value |
| --- | --- |
| Source path | `/atctools/Synopsys/tcad/T-2022.03/tcad/T-2022.03-SP2/Applications_Library/GettingStarted/swb/SimpleMOS/sdevice_des.cmd` |
| SHA-256 | `3fd501cad368aa43c55a780253309428872a4b3349aa03727f8a875d04a4b536` |
| SDevice banner | `Version T-2022.03-SP2` |
| Build identifier | `0.7745337, x86_64, Linux` |
| Verification date | `2026-08-26` |

The source deck itself and proprietary binary results are not added merely by
this evidence freeze. Generated TDR, PLT, log, and Vela state files remain
under ignored `build*/reference_tcad/simplemos_sentaurus2022` directories.

## Scope gates

- The input TDR is immutable and is not evidence for Vela process simulation.
- Carrier statistics are Maxwell-Boltzmann because the deck does not enable
  `Fermi`.
- OldSlotboom must run with the Fermi correction disabled.
- `SRH(DopingDependence)` is a recombination lifetime selection, not a
  low-field mobility selection.
- `Mobility(PhuMob HighFieldSaturation Enormal)` must retain its Mobility
  scope in the execution IR.
- Until PhuMob is implemented, production import must fail closed. It must not
  silently substitute Masetti or another low-field mobility model.
- Sentaurus normalized quasistationary step controls must not be copied into
  Vela as physical voltage steps without conversion.
- The comparison coordinate is `gate OuterVoltage`; the response is drain
  `TotalCurrent`, normalized to A/um.

## Milestones

| Milestone | Deliverable | Gate |
| --- | --- | --- |
| M0 | This plan, reference README, and machine-readable validation contract | Source provenance and SDevice-only scope are frozen |
| M1 | Scope-aware Physics parsing and fail-closed PhuMob classification | SRH doping dependence cannot enable Masetti mobility |
| M2 | One nominal TDR import gate | Four contacts, valid material regions, complete doping, qualified units |
| M3 | Restarted equilibrium/drain/gate workflow | Each stage consumes only the previous accepted state |
| M4 | Controlled A0-A3 mobility comparisons | Exact bias lattice and predeclared comparison contract |
| M5 | T-2022.03 PhuMob scalar implementation | Formula, parameter, and unit-system tests pass |
| M6 | Self-consistent PhuMob and Jacobian coupling | Directional derivative checks pass |
| M7 | T-2022.03 Enormal closure | Spatial Enormal and mobility comparisons pass |
| M8 | Original-deck B1 comparison | Nominal two-curve gate, then 8 TDR x 2 Vd matrix |
| M9 | Regression and documentation freeze | Focused and full test suites pass |

## M0/M1 implementation boundary

This first implementation stops after M1. It does not create a lossy PhuMob
substitution, import a process-generated TDR, change a solver default, or set
post-hoc scientific acceptance thresholds.

## M2 completion evidence

The nominal Workbench process node 17 was generated only to obtain the
immutable upstream SDevice input. Its parameters are `Lg=0.25 um`,
`NWell=1e17 cm^-3`, `GOxTime=10 min`, and `LDD_Dose=1e14 cm^-2`.

The resulting `n17_fps.tdr` has SHA-256
`a39f6eb1ed1745fad892b4ee74524dfce6f75a9287baeb4383faebb11856ee00`.
The M2 gate passes only when the TDR has the exact four contacts, valid
Silicon/Oxide/Nitride material regions, explicit qualified coordinate and
doping units, and complete donor/acceptor coverage on every Silicon region.
The TDR binary and generated neutral exports remain ignored build artifacts.
