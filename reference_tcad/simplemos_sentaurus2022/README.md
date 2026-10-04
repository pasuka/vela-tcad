# SimpleMOS Sentaurus 2022 SDevice reference

This directory defines the checked-in, SDevice-only comparison contract for
the Sentaurus Applications Library `GettingStarted/swb/SimpleMOS` project.

Sentaurus Process is outside the validation scope. Each process-generated TDR
will be treated as an immutable device input. SVisual plotting and extraction
are also outside scope; the comparison uses gate outer voltage and drain total
current directly.

## Current engineering entry (2026-09-30)

The [2026-10-01 raw asset archive and main migration](raw_assets_20261001/README.md)
adds persistent VM storage for local H5/TDR, solver logs and existing bundles,
with complete member hashes and source-to-main commit mapping.

Historical diagnostic files are indexed in the
[evidence archive and recovery guide](local_evidence_20260930/README.md).
Small records are checked in; detailed tables and figures have a fully hashed
persistent archive on the Sentaurus VM, independent of this disposable worktree.

Use [engineering/README.md](engineering/README.md) for the current validated
configuration, numerical run commands, external input requirements, and
byte-verifiable small evidence package. The n17–n24, Vd=0.05/1 V,
Vg=0:0.05:2.5 V matrix passed **816/816 dual-initialization points** with a
maximum absolute M60 current error of **0.186617394672%** against the 2% gate.
The subsequent Release regression passed **993/993**. The explicit profile
includes PhuMob, Enormal, HFS, signed Si source volumes and complete split state.

```text
python reference_tcad/simplemos_sentaurus2022/engineering/verify.py
```

This command only audits the shipped evidence. New simulations require the
external frozen meshes, complete split seeds and M60 native exports described
in the engineering README. No proprietary TDR/PLT or generated states are shipped.
The earlier M0–M4 contracts and commands below document historical stages;
their fail-closed model statements and earlier tolerances are not the current
engineering run instructions. Original inputs and historical evidence remain
unchanged. See also the [current branch status](../../docs/validation/simplemos_branch_status.md).

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

## Historical milestone status (M0–M4)

M0 freezes provenance and scope. M1 adds scope-aware SDevice model parsing and
keeps PhuMob fail-closed until an exact solver path is available. M2
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

M5 implements the T-2022.03 PhuMob scalar equations 270--283 and their Silicon
parameter set. The kernel accepts separate donor, acceptor, electron, and hole
concentrations plus temperature, and reports lattice, nonlattice, screening,
and combined mobility components. Golden formula tests and SI/TCAD-unit
invariance tests are frozen by
`simplemos_m5_phumob_scalar_contract_v1.json` and recorded in
`simplemos_m5_phumob_scalar_evidence.json`.

M5 deliberately does not make `Mobility.PhuMob` importable. Parameter coverage
uses the explicit status `scalar_kernel_only`, and production import remains
fail-closed until M6 supplies DD assembly and Jacobian coupling.

M6 connects the scalar kernel to both Gummel and coupled-Newton DD transport.
Each transport edge retains arithmetic endpoint averages of ionized donors,
ionized acceptors, electrons, and holes; compensated impurities are never
collapsed to net doping on the PhuMob path. Contact-current and fixed-state
diagnostics use the same state contract as the production residual.

The coupled Newton Jacobian differentiates the mobility response with respect
to both endpoint electrostatic potentials and both carriers' endpoint
quasi-Fermi potentials. This includes the cross-carrier blocks introduced by
carrier-carrier scattering. `Mobility.PhuMob` and
`Mobility(PhuMob HighFieldSaturation)` are now production-importable as
`phumob` and `phumob_field`. The original deck's additional `Enormal` selection
remains fail-closed until M7 supplies the T-2022.03 spatial closure; the older
Masetti/Lombardi control model is not substituted for it.

The machine-readable M6 contract and verification record are
`simplemos_m6_phumob_dd_contract_v1.json` and
`simplemos_m6_phumob_dd_evidence.json`.

M7 completes the T-2022.03 `EnormalDependence` path. The official Silicon
parameter export is frozen by SHA-256, and all 16 electron/hole parameter rows
are mapped in their native cm-based formula units before conversion to Vela's
active unit system. `Enormal` is composed with PhuMob as a low-field inverse-
mobility contribution, then the result is passed to high-field saturation.

The spatial closure selects only semiconductor/insulator interfaces, caches
interface normal and distance, and reconstructs the projected normal electric
field from the live cell potential. The coupled Newton pattern and Jacobian
therefore include the off-edge third vertex of each adjacent Tri3 transport
cell. Gummel assembly, contact current, and fixed-state audit use the same
geometry and state path.

The production importer now maps the original
`Mobility(PhuMob HighFieldSaturation Enormal)` selection to
`phumob_field_lombardi` with quasi-Fermi-gradient high-field drive. The M7
analytical spatial field, composed-mobility, and full-Jacobian comparisons all
pass. This makes the original deck executable for M8; it does not claim its
Id-Vg current parity yet.

The machine-readable M7 contract and verification record are
`simplemos_m7_enormal_contract_v1.json` and
`simplemos_m7_enormal_evidence.json`.
