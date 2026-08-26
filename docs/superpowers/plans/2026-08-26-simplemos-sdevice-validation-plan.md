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

## M3 completion evidence

`scripts/run_simplemos_m3_workflow.py` materializes two independent SDevice
branches at `Vd=0.05 V` and `Vd=1.0 V`. Each branch is the strict sequence
equilibrium -> drain ramp -> gate sweep. The equilibrium stage uses Vela's
Poisson-block initialization, followed by a coupled DD equilibrium point.

The workflow converts Sentaurus normalized sweep-time controls to Vela
physical-voltage controls by multiplying `InitialStep`, `MinStep`, and
`MaxStep` by the stage voltage span. It does not manufacture an exact gate
bias lattice; that remains an M4 deliverable.

Execution is fail-closed. A stage is accepted only when the runner succeeds,
all output rows are converged, the terminal voltage is reached, and a final
state file exists. The accepted state's SHA-256 must still match immediately
before its sole successor starts. This gate validates orchestration and state
continuity, not original-deck PhuMob physics or current parity.

## M4 controlled-comparison evidence

M4 defines the cumulative A0--A3 ladder `constant`, `masetti`,
`masetti_field`, and `masetti_field_lombardi`. The corresponding Sentaurus
decks use no explicit mobility selection, then add `DopingDependence`,
`HighFieldSaturation`, and `Enormal` one at a time. Carrier statistics,
OldSlotboom, SRH doping dependence, temperature, electrodes, TDR hash, and
solve sequence remain fixed.

The acceptance coordinate is exactly 51 direct gate points from 0 V through
2.5 V in 0.05 V increments at both drain voltages. The comparator rejects
missing, duplicate, non-finite, non-converged, or interpolated points before
evaluating its predeclared log-ratio, relative-error, endpoint, and trend
metrics.

The Vela continuation protocol keeps SRH density coupling, quasi-Fermi
reference, and continuity-row scaling invariant throughout each restarted
state chain. The default 1 mV auxiliary drain path, plus the bounded 5 mV
A3/1.0 V recovery described below, is allowed solely to reach the drain
operating points and is explicitly excluded from acceptance. The comparison
lattice and numeric thresholds were not changed after observing results.

The final qualification produced eight qualified Sentaurus reference curves
and eight Vela candidate curves, with 51 direct points per curve. All eight
comparisons pass. The matrix-wide worst absolute log10 ratio is 0.0455 dex,
the worst relative error is 11.1%, the largest absolute endpoint log10 ratio
is 0.0246 dex, and every curve has a matching nondecreasing trend.

The first A3/1.0 V pilot process exited with code `0xFFFFFFFF` after 302
accepted 1 mV auxiliary drain points and no Newton failure. Its bounded
recovery used 201 points at 5 mV, reached 1.0 V, and then converged all 51 gate
points. This recovery changes neither an acceptance point nor a comparison
threshold. M4 is complete; PhuMob remains excluded and fail-closed for M5.
