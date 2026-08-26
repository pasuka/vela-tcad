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
is recorded in `simplemos_nominal_tdr_m2_evidence.json`. No claim of SimpleMOS
current parity is made at M2.
