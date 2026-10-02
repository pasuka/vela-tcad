# WP5 orchestration validation, 2026-10-02

WP5 and WP6.1 are complete in `sdevice-cmd-frontend` on
`codex/sdevice-cmd-frontend`. WP6.1 is the separate commit `323983e0`;
its metadata and byte-compatibility evidence is in
[the FieldData record](sdevice_field_data_validation_2026-10-02.md).
The previous local-main merge and WP4 remain `cf10b3f1` and `5ae3ce13`.

## Implementation

- The `device` CLI resolves literal / `.cmd` / `_des.cmd` names and provides
  the frozen C8 options, numbered failures and exit codes 0 through 8.
- C1 parser, file existence, C2 mesh/electrode qualification, C3 materials,
  C4 decks and C10 sidecar gates run before any solver. `--check/-n` writes a
  check-only manifest with zero runner invocations; `-i` runs equilibrium only.
- Exclusive locks record PID, host, time and ownership token. Live/foreign or
  unreadable owners are preserved; dead local owners are reclaimed with W701.
  SIGINT/SIGTERM/Windows CTRL_BREAK stop the child, retain partial outputs,
  remove the owned lock and return 8. Windows liveness queries use process
  handles rather than `os.kill(pid, 0)`.
- Sequential stages use HDF5 state chaining. Partial stages without
  ExitOnFailure rebase later decks from the accepted CSV voltage. Save/Load
  snapshots carry measured partial voltages. Current Time points are bound
  from the preceding accepted current before execution; pure WP4 callers
  retain their previous default behavior.
- Complete, converged CSV rows are copied to a stable publication snapshot;
  PLT and PVD use atomic replacement. Finished VTU files are copied according
  to C4 output controls without editing binary payloads. The final file comes
  from the last accepted point. Quiet, final-only, binary and intermediate
  cleanup options are covered.
- ASCII runner inputs and temporary runtime log names avoid the runner's
  narrow Windows path conversion; Python handles Unicode paths and merges
  runtime records into the single user log at stage completion. CLI and logs
  use UTF-8. Manifests record decisions, overrides, stages and outputs;
  failures acquired under the lock also record `fail_report.json`.
- The CMake `vela_sdevice_pyz` target packs the package, four repository
  modules, schemas, default materials and current Git revision into
  `<build>/bin/vela_sdevice.pyz`; `vela` depends on it. The archive is usable
  without the source checkout, PYTHONPATH or third-party packages.

## Verification

MSYS2 UCRT64 GCC/C++20, Python 3.14, `windows-ucrt64-release`; HDF5 state
support and native TDR import were detected and enabled. Configuration also
detected the optional direct-solver packages; the imported test decks retain
the existing solver/model defaults. No Python API build was needed.

- Release configure and full build passed; native launcher and archive built.
- Metadata/VTU focused CTest: **4/4 passed**, 1.63 s.
- Full `ctest --preset windows-ucrt64-release -j 8 --output-on-failure`:
  **1020/1020 passed**, 156.04 s. Includes numerical and template regressions.
- SDevice CTest after final packaging/progress changes: **10/10 passed**,
  27.30 s. The orchestrator suite contains **19 tests**, with Windows signal,
  native launcher and real runner enabled.
- After persisting the electrode qualification findings, orchestrator,
  contracts and parser CTest passed **3/3**, 24.03 s (all 19 orchestrator cases).
- First full run was 1019/1020: a deliberate output-permission error exposed
  Windows' localized stderr encoding in the CLI capture. UTF-8 output fixed
  it; the complete rerun above passed. Physics and tolerance values were not
  changed to address this failure.

The fake runner verifies gates and all exit-code classes, conflicting/stale
locks, partial continuation, ExitOnFailure, abnormal/rejected children,
writer failure, signal cleanup, accepted-point progress and binary/final-only
publication. Hermetic tests clear PYTHONPATH and use isolated Python for
version, defaults, check-only gates and a complete fake-runner pipeline, with
site packages disabled (`python -I -S`). Explicit Plot/Time filenames are
reserved so the default accepted-point sequence cannot overwrite them.

Real-runner smoke uses the checked-in PN2D IV cmd fixture with only Grid,
Parameter and contact names adapted to a synthetic DF-ISE mesh: four nodes,
two triangles, 1 um square, uniform donor density `1e15 cm^-3`, no acceptors,
300 K. The fixture's model/Math/voltage/Plot controls are retained and
`--allow-lossy` is explicit. It creates log, PLT, VTU and PVD inside a Windows
directory containing Chinese characters and a space. A second real-runner
smoke uses Newton for equilibrium / voltage / Current-Time / HDF5 Load,
verifies all four stages reach their goals and restores the saved 0 V bias.
These are protocol/lifecycle checks, not PN2D numerical agreement evidence.

## Remaining boundaries

Current sign remains +1 and multi-stage PLT time resets per stage, explicitly
logged as pending WP10. GoalReached is local to each stage; the run manifest
records overall status. Linux signal/path smoke and distribution are WP11;
Sentaurus numerical comparisons and VM reference collection remain WP10/WP14.
No VM run, numerical calibration or generated simulation output is committed.
CLI errors before locating a cmd and lock conflicts do not overwrite another
run's files or create a misleading owned-workdir failure report.
