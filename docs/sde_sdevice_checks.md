# Generated SDE grids through sd-sim checks (P2)

P2 binds the three P1 mesh cases to eight named SDevice inputs, using the real
`sd-sim --check` launcher. It validates import, parameter/material identities,
electrodes, stage emission and sidecars. It does **not** run a numerical solver,
requalify an IV/BV curve, fetch private reference data or convert through TDR.

The frozen P2 contact offsets use the CODATA2018 constants profile. This
batch explicitly passes `--physical-constants-profile codata2018`, records that
profile in its summary and checks the emitted material identity. It does not
change the frontend default (`sentaurus_u2022_12`) or requalify those historical
offsets under a different constants profile.

The reviewed inputs and portable LF fingerprints live in
[`configs/reference_sde_checks.json`](../configs/reference_sde_checks.json), with
[schema](../schemas/vela.reference_sde_checks.v1.schema.json). The producer-side
geometry/doping contracts remain in the [P1 corpus](reference_sde_cases.md).

## Run from the worktree root

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --target sd-sim
python scripts/sentaurus_sde_cases.py `
  --output-dir "build-release/P1 meshes" --runtime-dir build-release/gmsh-python-runtime
python scripts/check_sde_sdevice_cases.py `
  --meshes-dir "build-release/P1 meshes" `
  --output-dir "build-release/P2 fresh checks" --launcher build-release/sd-sim.exe
```

For Linux, use its separate build directory, selected Python and built `sd-sim`
path; the same two Python CLIs apply. P2 requires a fresh empty output directory
and refuses to replace an earlier run or user files. Keep the existing directory
as evidence and choose a new one on the next invocation.

Exit 0 plus `check_summary.json.status=passed` means all eight check-only inputs
met the contracts. Exit 2 from this batch tool means an identity or check failed;
individual native command exit codes and their original logs remain in the
local case directories. The wrapper does not change native exit-code semantics.

## Named inputs and explicit staging

| Mesh | SDevice variants | Parameter/sidecar policy |
| --- | --- | --- |
| PN2D | 0V, IV, BV | Copy the repository PN2D `models.par`; BV also copies the exact existing atomic-profile sidecar |
| gate_reference | Barrier zero/positive/negative, WorkFunction | Copy its declared reference `models.par`; preserve all four command boundary values |
| Schottky | `schottky_n_des.cmd` | Explicitly stage the current Vela default `models.par` **for check-only legality**, not as a claim of native-reference parameter parity |

Each command is copied into an isolated case directory. Only `File{Grid}` and
`File{Parameter}` declarations are rebound to `device.grd` and `models.par`.
The missing Parameter declaration is explicitly inserted in the copies that need
it. Original files remain unchanged; their raw/LF hashes, effective-command hash,
edits, original parser verdict and material-file policy are recorded. In
particular, an original PN2D command with E301 is not relabelled as accepted.
No Physics, electrodes, goals or frozen input file is edited.

The tool compares the original and effective parsed electrodes, physics, plots,
Math, Solve and stage statements. It binds each P1 generated-file hash before
launching the frontend, and rechecks the analytic P1 geometry/doping/quality
contract. Missing or modified P1 inputs fail before native invocation.

## Interface import correction

DF-ISE `material=Interface` records are **line boundaries**, not volume regions
or material-library entries. The importer preserves their name, nodes, complete
edge segments and adjacent volume-region names as optional `mesh.json.interfaces`
metadata. Every interface edge must separate two volume regions; non-line and
unattached records fail E402. Solver regions and the C2 `materials` list contain
only actual volume materials. The C2 qualification schema and numerical policies
are unchanged; the C++ mesh reader does not use this optional boundary provenance.

This fixes the observed gate-reference E403 caused by requesting a fictitious
`Interface` material from the parameter library. Gate contacts remain attached to
the oxide; they are not converted into silicon contacts or ohmic electrodes.
P1 identity comparison now also checks the preserved interface edge set.

## Verified identities and boundaries

- All eight inputs must pass the native check with `check_only=true`, no executed
  stages, no state files or stage CSVs. A solver sentinel detects an accidental
  numerical invocation; the imported source must be DF-ISE with matching GRD/DAT
  hashes, so no TDR conversion is part of this path.
- The copied Parameter hash, runtime default material-library hash, generated
  material hash, launcher/zipapp fingerprints and source/effective command hashes
  are bound. There is no automatic `--allow-lossy` escape.
- Coordinates stay in um; doping in cm^-3 and energies in eV. Region-labeled
  topology, positive orientation, contact nodes/edges/owners, interface edges and
  donor/acceptor/net values are compared to the producer bundle.
- Emitted stage kinds, sequential dependencies, state filenames, all electrode
  start/end voltages and actual JSON sweep goals must match the source statements.
  Merely producing a valid `stage_plan.json` is insufficient.
- Schottky thermionic contact type, work function and electron/hole surface
  velocities are preserved; gate Barrier/WorkFunction mapping remains the current
  qualified frontend mapping. PN2D BV retains the complete SG/GSS-Laux source,
  Bernoulli, mixed-Voronoi and non-obtuse bundle plus the declared field recovery.
- Actual negative controls retain E301/E303 for Parameter failures, E401 for an
  unmatched electrode, E403 for an unknown volume material, E106 for unsupported
  Thermodynamic physics and E409 for a sidecar trying to replace Physics.

## Qualification and tests

```powershell
ctest --test-dir build-release --output-on-failure -R '^sde_sdevice_checks$'
```

CTest passes its actual `sd-sim` target through `VELA_SD_SIM`. Direct unittest
runs may set that variable explicitly; developer machines without Gmsh or the
launcher skip the native integration class visibly. The batch CLI fails if its
required launcher or grids are missing.

The current numerical qualification of PN2D remains tied to the named
1943-node/3680-cell reference grid and its inputs/materials/profile in the
[delivery matrix](sdevice_support_matrix.json). P2's generated 1482/2800 grid is
a separate experimental arm. This work retains the frozen control identity in
its acceptance record and does not recompute or supersede its numerical evidence.
Gate-reference uses the user-approved P1 contract-equivalent BAMG policy; different
mesh hashes remain visible. Full numerical controls, complete IV/BV and Linux
native qualification are separate gates.

Keep GRD/DAT/MSH, logs, TDR and command copies in ignored build directories.
Commit only the small contracts, source/tests/docs and structured acceptance
evidence. The next planned implementation is P3: BJT body handles, General
analytical doping and bounded MaxTransDiff refinement.
