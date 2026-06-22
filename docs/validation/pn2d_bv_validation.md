# PN2D BV Validation Methodology

## Source Of Truth

The source deck is
`reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd`. It uses
`Avalanche(VanOverstraeten)` and sweeps Anode to `-20.0 V`.

## Artifact Refresh

Sentaurus artifacts are refreshed through
`scripts/run_sentaurus_vm_reference.py pn2d --stages 0v,iv,bv`. Dry-run planning
uses `--dry-run`; live upload/run/fetch omits `--dry-run`.

## Required Artifact Checks

The BV refresh must contain `pn2d_bv.plt`, a clean BV log artifact such as
`pn2d_bv.log` or Sentaurus-generated `pn2d_bv.log_des.log`, and 201
endpoint-inclusive `pn2d_bv_multibias_*_des.tdr` files.

## Model Decision

When the refreshed source remains `Avalanche(VanOverstraeten)`, Vela's existing
`van_overstraeten` implementation is the target and Okuto-Crowell remains
contrast-only.

The refreshed Sentaurus BV source and currently available generated artifacts
use `Avalanche(VanOverstraeten)`. No Okuto-Crowell model is added in this
validation pass.

## Comparison Layers

1. Curve comparison checks current trend and documented windows.
2. Field comparison checks potential, electric field, carrier density, mobility,
   and avalanche generation/source density at selected biases.
3. Coefficient checks compare alpha to alpha; generation checks compare
   generation/source integral to generation/source integral.

## Accepted Gates

The promoted automated gates are VM-free and lightweight: parser/import
provenance, artifact validation with synthetic fixtures, BV max-field trend, and
documented comparison summaries.

## Non-Goals

This pass does not claim full `0..-20 V` absolute-current parity, does not
promote hidden scalar source calibration, does not rewrite SG flux divergence,
and does not add LDMOS/IGBT/MOS BV validation.

## High-Bias Interpretation

Windowed current agreement and high-bias knee shape are reported separately. If
knee-shape evidence remains divergent, the methodology records it as an open
physics/parity limit rather than hiding it behind a broad current band.

## Current Diagnostic Result

The derived Vela `-20.0 V` candidate deck under
`build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_validation_candidate`
currently records a controlled non-convergence diagnostic before the requested
high-bias comparison points. The latest local run reached a last stable bias of
about `-0.2301056002 V` and then failed at `-0.2301056003 V` with
`max_iterations`. The generated comparison report therefore treats the requested
`-0.5, -2, -5, -10, -20 V` candidate curve and field points as missing and
leaves the `-13.2..-13.0 V` current window `not_evaluated`; this is diagnostic
evidence, not a promoted `-20 V` parity gate.

After the live `pn2d_bv_validation_refresh` Sentaurus VM run, the refreshed BV
artifacts validate with 201 endpoint-inclusive multibias TDR files and a clean
BV log artifact. The local Vela candidate still does not reach `-20.0 V`.
Manual cross-checks against the Sentaurus Device user guide led to one code
alignment in this pass: GradQuasiFermi avalanche coefficient fields now fall
back to electric field in contact-touching elements. Control runs show that this
is not sufficient for the high-bias continuation limit: `electric_field` driving,
legacy `mobility_density_gradient`, `impact_ionization = none`, and a higher
Newton iteration budget all still stop around `-0.20..-0.23 V`. The remaining
gap is therefore treated as a coupled DD continuation/contact-state problem, not
as resolved BV physics parity.
