# Changelog

This log starts with the SDevice frontend delivery. Earlier project changes
remain in Git history and the documentation index.

## Unreleased — SDevice frontend closeout (2026-10-05)

- Aligned built-in Si/SiO2 scalars to the emitted no-BGN 300 K material library,
  including self-consistent intrinsic density. Added a generated snapshot,
  dual-unit C++ checks and explicit legacy materials for calibrated reference
  workflows; original partial inputs and numerical tolerances are preserved.

- Added `vela device` for the documented preprocessed single-device SDevice
  subset: input gates, mesh/material conversion, v2 stage decks, HDF5 state
  chaining, `--check`, default parameter export and numbered failures.
- Added accepted-point PLT/VTU/PVD publication, per-stage Time and goal metadata,
  partial outputs on interruption, owned locks and controlled solver threads.
- Added offline relocatable Windows/UCRT64 and Ubuntu 26.04 x86_64 packages with
  CPython 3.12, runtime libraries, full license notices and file/commit manifests.
  Retired the earlier Ubuntu 24.04 environment and artifacts.
- Repaired PN2D IV/BV model activation, contacts and continuation. Windows native
  reference acceptance reaches 10/-20 V; curve and spatial qualification remain
  limited to the explicitly pinned contract. The -13.2 V spatial sample is not
  the breakdown voltage; near-19.6 V diagnostics and the named curve-knee metric
  are documented, and the -20 V hole-density difference remains unresolved.
- Completed workflow/schema/import documentation, supported-subset and message
  audits, executable neutral mesh/sidecar examples and the six-item DoD ledger.
  Fresh native PN2D 0V/IV/BV acceptance and check-only gates pass through both
  extracted packages under the pinned contract, with unchanged package file
  fingerprints. WP13 built-in alignment is recorded separately; the remaining
  native reruns (WP14) are independent follow-ups.
- Fixed stage sidecar `runtime_log.profile` being overwritten by the
  orchestrator. CLI verbose selects debug globally; otherwise each stage uses
  its sidecar profile or default. Real runner log regression covers precedence.
- Fixed Windows license collection for CMake preset compiler cache entries and
  disabled bytecode writes explicitly in the isolated native E2E controller.
  The invalid-cache regression passes on both platforms; original failed
  package integrity evidence is retained in the follow-up record.

See [workflow](docs/sentaurus_sdevice_workflow.md),
[packaging](docs/sdevice_packaging.md) and
[WP12 evidence](docs/validation/sdevice_wp12_closeout_2026-10-05.md) and the
[native package follow-up](docs/validation/sdevice_native_package_e2e_2026-10-05.md).
