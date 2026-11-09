# Sentaurus SDE/SDevice end-to-end acceptance

The public entry point is:

```powershell
python scripts/sentaurus_import.py device `
  --sde <sde.cmd> --sdevice <sdevice.cmd> `
  --output-dir <output> --sde-frontend offline-static --run `
  --runner <absolute-path-to-build-release/vela_example_runner.exe>
```

`offline-static` is also the CLI default, so the explicit option above is
documentary rather than required. Public distributions and CI must remain
license-free and must not invoke Sentaurus. The
`sentaurus2022-evaluated` frontend is restricted to internal development and
reviewed oracle generation on licensed infrastructure; it is never a runtime
dependency for end users.

The command emits frontend provenance, Device IR, Execution IR, mesh, doping,
mesh qualification, Vela deck, curve, acceptance report, and run manifest.
Inputs matching a case in `configs/sentaurus_e2e_contract.json` are evaluated
against that case automatically; `--contract-case` selects one explicitly.

`generated` means artifacts and deck exist but no solve ran. `solved` means the
runner completed but no unique contract case matched. `accepted` means the
runner and every contract gate passed. `failed` means parsing, generation,
solve, or acceptance failed; `fail_report.json` identifies the stage and, for
curve coverage failures, the last accepted bias.

The full PN2D BV case is deliberately not declared complete merely because a
partial curve exists. If it cannot reach the contracted stop voltage, the
acceptance gate fails with the actual last bias.

For the longer full-mesh PN2D cases, use a Release runner from an independent
`build-release` directory. Do not reuse Debug object files in that directory.
The Python import and acceptance tests remain valid against either runner.
