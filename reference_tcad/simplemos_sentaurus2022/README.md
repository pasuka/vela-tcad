# SimpleMOS Sentaurus 2022 SDevice reference

This directory defines the checked-in, SDevice-only comparison contract for
the Sentaurus Applications Library `GettingStarted/swb/SimpleMOS` project.

Sentaurus Process is outside the validation scope. Each process-generated TDR
will be treated as an immutable device input. SVisual plotting and extraction
are also outside scope; the comparison uses gate outer voltage and drain total
current directly.

## Frozen source

- release: `T-2022.03-SP2` (`0.7745337, x86_64, Linux`)
- source: `/atctools/Synopsys/tcad/T-2022.03/tcad/T-2022.03-SP2/Applications_Library/GettingStarted/swb/SimpleMOS/sdevice_des.cmd`
- source SHA-256: `3fd501cad368aa43c55a780253309428872a4b3349aa03727f8a875d04a4b536`
- live verification date: `2026-08-26`

The original deck selects:

- four electrodes: source, drain, gate, and substrate;
- Maxwell-Boltzmann carrier statistics;
- `EffectiveIntrinsicDensity(OldSlotboom)`;
- silicon `Mobility(PhuMob HighFieldSaturation Enormal)`;
- silicon `Recombination(SRH(DopingDependence))`;
- a Poisson initialization, coupled DD equilibrium, drain ramp, and gate
  sweep to 2.5 V.

## Repository policy

Only neutral text inputs, contracts, hashes, configurations, comparison
tables, and reports may be checked in. Proprietary TDR/PLT files, Sentaurus
logs, generated Vela states, and plots remain below ignored build trees.

The machine-readable contract is
`simplemos_sdevice_validation_contract_v1.json`. Generated evidence should be
staged below `build-release/reference_tcad/simplemos_sentaurus2022`.

Run the M2 gate with:

```text
python scripts/run_simplemos_m2_tdr_gate.py --tdr <n17_fps.tdr> --output-dir <ignored-output-dir>
```

The command verifies the frozen SHA-256 before invoking `sentaurus_import`, so
a different binary cannot pass merely because it has a compatible topology.

Materialize the M3 SDevice workflow with:

```text
python scripts/run_simplemos_m3_workflow.py --base-config <qualified-vela-device-deck.json> --output-dir <ignored-output-dir> --clean
```

The default branches use `Vd=0.05 V` and `Vd=1.0 V`. Each branch owns an
independent equilibrium, drain-ramp, and gate-sweep state chain. Add
`--execute --runner <vela_example_runner>` to run a physics-ready base deck.
Workflow acceptance proves restart continuity only; it does not waive the
fail-closed PhuMob gate or claim current parity with the original deck.

## Current status

M0 freezes provenance and scope. M1 adds scope-aware SDevice model parsing and
keeps PhuMob fail-closed until an exact implementation is available. M2
qualifies the nominal node-17 process TDR as an immutable SDevice input:

- SHA-256 `a39f6eb1ed1745fad892b4ee74524dfce6f75a9287baeb4383faebb11856ee00`;
- 1,542 vertices, four material regions, four contacts, and 132 datasets;
- explicit `cm` geometry converted to Vela's neutral `um` coordinate export;
- complete scalar `cm^-3` donor and acceptor coverage on the Silicon region;
- zero exported doping on nodes that belong only to Oxide or Nitride.

The binary TDR remains in an ignored build directory. The neutral M2 evidence
is recorded in `simplemos_nominal_tdr_m2_evidence.json`.

M3 adds a machine-checkable restarted workflow for both target drain voltages.
Sentaurus sweep-time controls are scaled by the physical voltage span: the
0.05 V drain branch starts at 0.005 V, the 1.0 V branch at 0.1 V, and both gate
sweeps start at 0.025 V with a 0.125 V maximum step. A successor is launched
only after the predecessor reaches its terminal bias with every emitted point
converged and its final state hash recorded. No claim of SimpleMOS current
parity is made at M3.

The neutral M3 implementation and regression hashes are recorded in
`simplemos_m3_workflow_evidence.json`; generated manifests and restart states
remain ignored build artifacts.

M4 freezes a controlled mobility ladder that deliberately excludes PhuMob:

| Variant | Sentaurus mobility | Vela mobility |
| --- | --- | --- |
| A0 | default constant mobility | `constant` |
| A1 | `DopingDependence` | `masetti` |
| A2 | A1 + `HighFieldSaturation` | `masetti_field` |
| A3 | A2 + `Enormal` | `masetti_field_lombardi` |

All eight reference curves use `Vd=0.05 V` or `1.0 V` and the exact direct
gate lattice `0:0.05:2.5 V` (51 points). Interpolation is forbidden. The Vela
workflow uses extra drain-ramp points only for continuation; these points have
no comparison or acceptance role. Solver settings that affect the state-chain
equations remain invariant across equilibrium, drain, and gate stages.

Prepare or execute the matrix with:

```text
python scripts/run_simplemos_m4_controlled_matrix.py --execute-vela
```

Run the fail-closed exact-grid comparison with:

```text
python scripts/compare_simplemos_m4_controlled_matrix.py --output-dir <ignored-report-dir>
```

The checked-in `controlled_mobility` directory contains only the qualified
neutral 51-point reference CSVs and their manifest. Proprietary raw files and
Vela run outputs remain below ignored build directories.

Final M4 qualification passes all eight curves and all 408 direct points.
Across the matrix, the worst absolute log10 current ratio is 0.0455 dex, the
worst relative error is 11.1%, and all eight current-magnitude trends match.
These results pass the frozen 0.3 dex, 100% relative-error, and trend gates.
The complete per-case metrics and the bounded A3/1.0 V recovery record are in
`simplemos_m4_controlled_mobility_evidence.json`.

The device mesh, net-doping, Id-Vg, and direct nodal-field figures are
documented in `../../docs/validation/simplemos_m4_visual_report.md`. Their
reproducible comparison inputs are under `controlled_mobility/comparisons`.
