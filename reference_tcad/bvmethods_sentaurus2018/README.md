# BVMethods NMOS reference family

## Current qualified snapshot: 2026-11-09

The `bvmethods_sentaurus2018` directory name is historical. Use
[current_validation.json](current_validation.json) and the versioned
[U-2022.12 compatibility snapshot](u2022_12_20261109/README.md) for the current
qualified Vela BVMethods CMD/JSON defaults and evidence. See the
[full closure report](../../docs/validation/bvmethods_production_closure_2026-11-09.md).

| Dataset | Scope | Accepted / attempts / failed |
|---|---|---:|
| Native Sentaurus U-2022.12 | Adaptive 0-100 V only | 45 / 48 / 3 |
| Vela ordinary physical-coordinate entry | Full BV to current stop | 62 / unknown / unknown |
| Vela ordinary Unit-coordinate entry | Full BV to current stop | 62 / unknown / unknown |
| Vela original CMD stage chain | Full BV to current stop | 106 / unknown / unknown |

This snapshot uses 300 K, the 2600-node/5020-Tri3 mesh, matching U-2022.12
constants/materials, SRH plus GradQuasiFermi avalanche, BTBT off and a
1e7 ohm*um resistor. Candidate curves are under the version's `vela/curves/`;
the native bounded curve is under its `reference_curves/`. Current CMD/par,
sidecar, material data, relative stage configurations and input SHA receipts
are included; binary mesh and bulk generated mesh/doping remain external.

The native bounded curve does not qualify a native full BV. Vela accepts use
the explicitly approximate E-or-mapped-RHS rule; strict native stopping and
performance noninferiority remain unqualified. Carrier/continuity report
failures are retained, including 70/13 for the CMD chain. Different entry/stage
chains prevent treating its larger workload as a coordinate-only comparison.

This is a local text qualification snapshot. It does **not** replace the native
artifact-store `artifact_manifest.json` (still v4), its v1-v4 predecessors,
`source/`, `wp14_sources/`, legacy materials, reference curves or historical
reference inventories. Their existing SHA contracts and physics remain intact.
Dated investigations below are provenance, not active task queues.

## Historical O-2018 and October references


On 2026-10-15, the user-selected target for the external-resistor
investigation was the installed T-2022.03-SP2 library's original `resistor`
branch: NonLocalPath BTBT, GradQuasiFermi avalanche, and a `1.3e-3 A/um`
stop. See the [dated alignment plan](../../docs/superpowers/plans/2026-10-15-bvmethods-installed-library-alignment-plan.md).
The checked-in E2/Eparallel sources and existing qualifications retain their
historical scope; they have not been rewritten or requalified as library originals.
The [2026-10-15 original-reference audit](../../docs/validation/bvmethods_original_reference_q0_q1_2026-10-15.md)
reproduces the archived T-2022.03-SP2 curve byte for byte in a new from-zero run;
its native formatted BV is 6.67 V (independent linear interpolation: 6.673594819 V).
Vela still rejects NonLocalPath with E106; this native reference is not Vela equivalence.

This directory preserves the text inputs used by the Sentaurus Training
`BVmethods` NMOS example and the method mapping used by Vela validation. Binary
TDR/PLT/log outputs remain generated artifacts under `build*/reference_tcad/`
and are not checked in.

## Source inventory

| Checked-in file | Training workbench node | Method |
|---|---:|---|
| `source/bvmethods_nmos_sde.cmd` | 1 | symmetric NMOS structure and mesh |
| `source/aba_poisson_sdevice.cmd` | 3 | ABA, Poisson-only sweep |
| `source/aba_coupled_sdevice.cmd` | 4 | ABA postprocessing on coupled DD |
| `source/external_resistor_sdevice.cmd` | 5 | `1e7 ohm*um` series-resistor load line |
| `source/voltage_to_current_sdevice.cmd` | 6 | voltage ramp to 6 V, then current control |
| `source/continuation_sdevice.cmd` | 7 | Sentaurus continuation curve tracing |
| `source/transient_sdevice.cmd` | 8 | transient voltage ramp with resistor |
| `source/models.par` | 3--8 | shared SRH lifetime parameters |

The SDevice files preserve the training method mapping and contact names;
file/output names allow the decks to coexist outside Sentaurus Workbench.
They are a repository adaptation, not a verbatim copy of every archived node.
For example, archived `v1/source/pp5_des.cmd` declares
`Band2Band(Model=NonLocalPath)`, `Avalanche(GradQuasiFermi)` and a
`1.3e-3 A/um` stop, while the current external-resistor source declares
`Band2Band(E2)`, `Avalanche(Eparallel)` and `1.443e-3 A/um`.
The historical final current sample must not be confused with its declared
stop. See the [zero-QS diagnosis](../../docs/validation/bvmethods_zero_qs_diagnosis_2026-10-06.md)
for the historical E2/Eparallel contract and retained failure evidence.

Run the structure deck first. It builds `bvmethods_nmos_half_msh.tdr`, mirrors
it about `x=0`, renames the mirrored drain contact to source, and writes
`bvmethods_nmos_msh.tdr`. Each SDevice deck expects that final mesh in its
working directory.

## Vela template mapping

- `vela/materials_sentaurus2018.json` preserves the Si/Nitride material
  overrides used by the boundary-method and IIC replay scripts. The values
  are unchanged; these inputs are now stored with the BVmethods case.
- `configs/templates/bvmethods_nmos_external_resistor.template.json` maps the
  series-resistor load-line method.
- `configs/templates/bvmethods_nmos_voltage_to_current.template.json` maps the
  voltage-to-current boundary switch.
- IIC/ABA validation remains represented by the existing BV/IIC configuration
  workflow and the checked-in validation scripts.
- Sentaurus `Continuation` is mapped to Vela's existing pseudo-arclength BV
  mode. The five non-transient references and the current Vela closure status
  are sealed in `bvmethods_nontransient_validation_20260817.json`.
- Sentaurus `Transient` remains archived as a reference only and is explicitly
  outside the non-transient closure scope.

The v4 C7 inventory records SHA-256 and sizes of the exact archived source
bytes, including CRLF where present. Git preserves `source/` bytes accordingly;
these hashes differ from LF-normalised `producer.cmd_sha256`. The v4
`source/wp14_inputs.json` receipt is included with the authored inputs so its
repository provenance can be verified. Raw simulation results remain external.

The historical O-2018 reference threshold voltages are `6.379791636 V` for the
external-resistor method and `6.383184201 V` for voltage-to-current control at
the example drain-current criterion. No mobility or avalanche coefficient
scaling is used by the Vela templates.

`scripts/freeze_bvmethods_nontransient.py` regenerates the five compact
reference curves from the archived PLT files. The Continuation deck is
prepared by `scripts/prepare_bvmethods_nmos_continuation.py`, which copies an
accepted base deck's physics verbatim and changes only its sweep controller.
The preparer requires two adjacent accepted states: the earlier state supplies
the initial secant direction and the later state is the restart anchor. This
is a numerical initialization only; it does not modify the device equations.

The current NMOS continuation status remains **pending**. Analytic-tangent
trials with and without Enormal mobility advanced only about `1.5e-10 V`
before reaching the minimum arclength step. A bounded two-state secant trial
from 6.000000 V to 6.056459 V produced no accepted point after the anchor.
These results are sealed in `continuation_diagnostic_20260818.json`; they do
not block the already passed IIC, external-resistor, or voltage-to-current
branches and do not justify adding another physical model.
