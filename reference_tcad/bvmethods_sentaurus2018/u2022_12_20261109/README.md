# Qualified U-2022.12 compatibility snapshot (2026-11-09)

This text-only snapshot accompanies integration commit `77a889a9` and the
[closure report](../../../docs/validation/bvmethods_production_closure_2026-11-09.md).
The parent directory's 2018 name is historical. Its old sources, curves,
materials and artifact manifests remain frozen. Vela candidates are not
native references; this snapshot is not a new remote-store registration.

## Contents and qualification

- `source/` retains the exact validated CMD/par bytes and explicit basic-state
  Plot sidecar. Unsupported current/source Plot exports are not qualified.
- `vela/materials.json` and `simulation_stage00/01.json` are the actual matching
  material and portable stage configurations, without numerical retuning.
- `vela/curves/` contains ordinary JSON full-BV candidates (62 points each)
  and the original-CMD candidate (106 points). Successful sweep Newton counts
  are 351, 351 and 732. Initialization/stage chains differ.
- `reference_curves/` contains the native U-2022.12 adaptive 0-100 V curve:
  45 accepted / 48 attempts / 3 failures, 308 nonzero-target Newton iterations.
  The raw PLT has three initialization records at zero; the neutral curve
  retains the last record per outer bias, matching the qualified ledger.
  Native zero-bias iterations are unknown in this nonzero-goal accounting.
- `comparison.json` records scope, units, report deficits, identities,
  unknown counters, qualification limits and evidence hashes.

Candidates reach the 1.3e-3 A/um stop and pass applicable F/E, load-line,
KCL/finite checks and the existing BV interpolation gate. Carrier/continuity
report failures remain 3/7, 3/8 and 70/13. Performance noninferiority and strict
native RHSMin equivalence are unqualified. The BV anchor is an existing
acceptance contract, not a BV measured from the bounded native curve.

## Reproduction inputs

Binary mesh and generated bulk mesh/doping remain external. Verify the mesh
against `source/inputs.json`, placing `n1_msh.tdr` beside `pp5_des.cmd`.
Use a matching `sentaurus_u2022_12` build with HDF5 import and PARDISO/MKL.
Set `PYTHONPATH` to the absolute repository `scripts/` directory (or use an
installed frontend package). Select `VELA_LINEAR_SOLVER=pardiso_mkl`, set
`VELA_MKL_RT_LIBRARY` to your local runtime and
`MKL_THREADING_LAYER=SEQUENTIAL`. From this snapshot directory:

```text
python -m vela_sdevice device source/pp5_des.cmd --allow-lossy --threads 1
```

The frontend generates mesh/doping/materials and stages in `source/pp5_des_vela`.
For direct use of the saved `vela/` configurations, first stage the verified
`mesh.json` and `doping.csv` in that folder and change into `vela/`. Run `simulation_stage00.json`, then
`simulation_stage01.json`, which loads `stage00_state.h5`. Backend selection
uses the documented runtime configuration. Shared-VM wall time is not an
exclusive speed benchmark. Output paths are relative to the run directory.
Use an external 12-hour supervisor as in the qualification; these decks do
not install a wall-clock supervisor by themselves.

No numerical run was repeated for this synchronization. Regenerate with
`scripts/freeze_bvmethods_u2022_validation.py` from the repository root,
explicit validated evidence/native inputs and a fresh version output folder.
