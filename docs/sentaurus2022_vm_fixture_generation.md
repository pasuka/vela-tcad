# Generating authoritative Sentaurus 2022 evaluated fixtures

Run this procedure inside the licensed Sentaurus 2022 VM. The repository CI
does not require a license and consumes only reviewed, frozen artifacts.

1. Check out the same commit that will receive the fixture.
2. Put the Sentaurus T-2022.03 `sde` executable on `PATH`.
3. Run the paired workflow with the evaluated frontend, for example:

   ```powershell
   python scripts/sentaurus_import.py device `
     --sde reference_tcad/pn2d_sentaurus2018/source/pn2d_sde.cmd `
     --sdevice reference_tcad/pn2d_sentaurus2018/source/pn2d_iv_sdevice.cmd `
     --sde-frontend sentaurus2022-evaluated `
     --sde-build-mesh `
     --output-dir build/sentaurus2022-fixtures/pn2d-iv
   ```

4. Review `sde.stdout.txt`, `sde.stderr.txt`, exit code, SDE version and both
   hashes in `sde_frontend_manifest.json`.
5. Preserve the original owned `.cmd`, `.log.eval`, version/provenance manifest,
   expected `device_ir.json`, `_bnd.tdr`, `_msh.cmd`, and selected golden sample
   values. Do not commit license tokens, VM paths, unrelated generated plots, or
   private third-party decks.
6. Re-run the offline-static path against the frozen `.log.eval`; its Device IR
   must match the reviewed evaluated IR except for `source` provenance.

At minimum generate four fixture families: ABA/BAB multi-region overlap;
`insert-vertex` plus partial edge contact; constant Replace/LocalReplace/
NoReplace; and directed line-baseline Gaussian Positive/Negative/Both plus
NoEval. Hash changes require fixture review, not automatic golden replacement.

This development machine did not expose an `sde` executable on `PATH`, so no
official T-2022.03 `.log.eval` or TDR oracle has been fabricated locally.
