# PN2D Sentaurus 2022 SDE oracle

This frozen oracle was generated from
`reference_tcad/pn2d_sentaurus2018/source/pn2d_sde.cmd` by the licensed
Sentaurus Structure Editor T-2022.03-SP2 VM using
`scripts/run_sentaurus2022_sde_fixture_vm.py`.

`manifest.json` is the provenance root. It seals the owned input, evaluated
wrapper, official variable-expanded command log, boundary TDR, generated mesh
command, mesh TDR, diagnostics, version output and expected Vela Device IR.
The official mesh summary is 1943 points, 3680 elements and 3 regions.

The source is intentionally not modernized inside this fixture. The preserved
stdout contains T-2022.03's warning that `sdegeo:define-2d-contact` is obsolete
and should be replaced by `sdegeo:set-contact`; changing the owned source must
create a reviewed oracle revision with new hashes.
