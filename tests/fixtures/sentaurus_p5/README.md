# Owned P5 native controls

These short Scheme inputs were authored for Vela's P5 investigation. They do
not contain vendor implementation code or output data.

- `ratio_one.scm`: 0.32 um square, max=0.08/min=0.02 um, ratios=(1,1).
  Create two controls by setting min=max=0.08 and min=max=0.02 respectively;
  use distinct output basenames `ratio_max` and `ratio_min`.
- `body_remove.scm`: overlapping gate removal, a subsequent deleted-region
  doping reference, and a surviving silicon area of 0.82 um^2.
- `erf.scm`: a Gaussian normal profile with a finite lateral Erf baseline and
  an inserted endpoint vertex, including the half-amplitude endpoint.

Run them only in an existing licensed Sentaurus environment. Use the private
VM workflow in `docs/office_vm_workflow.md` for remote runs; do not publish host
or license settings. Convert outputs with `tdx -dd <mesh-basename> <output>`.
The auditor expects basenames `native_ratio`, `native_min`, `native_max`,
`native_removed`, `native_erf`, and `native_bv`. BV's unchanged primary source
needs a separate `+z` context wrapper as described in `docs/sde_bv_modeling.md`.

Keep every generated GRD/DAT/MSH/TDR and raw log ignored. The public structured
acceptance contains numeric statistics and hashes only.
