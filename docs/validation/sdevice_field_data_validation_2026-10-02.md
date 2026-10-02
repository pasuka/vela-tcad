# WP6.1 FieldData validation, 2026-10-02

Worktree: `sdevice-cmd-frontend`; toolchain: MSYS2 UCRT64 GCC, C++20,
`windows-ucrt64-release`. HDF5 state files and TDR import were detected and
enabled. This change only parses output metadata and supplies it to the
existing VTU writer; physics, discretization and numerical tolerances are unchanged.

`sweep.vtk_field_data` accepts ASCII identifier keys and string/finite number
values, rejects reserved keys and invalid value types, and enables per-point
time, local goal status and swept-contact voltages. Explicit custom `t` is
preserved; missing `t` is computed from actual stage voltage/current progress.
Both `GoalReached` and the C6 alias `goal_reached` are written only when the
object is present. Existing decks without it retain their previous output bytes.

The new Catch2 case independently reads ASCII and binary VTU metadata,
checks `CmdSHA256=abc`, `Stage=1`, explicit `t=0.25`, computed midpoint
`t=0.5`, initial/final goal status, invalid types/names/reserved keys, and
byte equality after removing the opt-in key. A four-node, two-triangle
Si device is used, with coordinates in um and donor density `1e15 cm^-3`.

Verification:

- Release configure and complete build: passed.
- `ctest --preset windows-ucrt64-release -R 'DCSweep frontend FieldData|VtuWriter'`:
  **4/4 passed**, 1.63 s; includes the new metadata case and independent VTU readers.
- `ctest --preset windows-ucrt64-release -j 8 --output-on-failure`:
  **1020/1020 passed**, 156.04 s, including the WP5 integration tests.
  The first unlimited-parallel run passed 1019/1020; WP5's intentional I/O
  failure exposed a Windows stderr encoding mismatch, fixed in its CLI before
  this complete rerun. No numerical regression failed.

No generated simulation output is committed. WP10 comparison with Sentaurus
and current-sign calibration are outside this validation.
