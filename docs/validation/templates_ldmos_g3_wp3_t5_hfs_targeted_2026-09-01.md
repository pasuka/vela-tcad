# Templates/LDMOS G3 WP3 T5 contract-equivalent H1 targeted control

Date: 2026-09-01
Frozen plan: `01f20ace88dc2b68a7bc3788534244dadf0d18f2` (R4)
Status: complete; main gate failed, migration guard passed, stop before T6

## Outcome

The pre-registered, single H1 candidate
`transport_cell_magnitude_average` was evaluated on the two frozen endpoints
by exact algebraic reconstruction of Vela's electron transport row.  It
preserves the complete Sentaurus state, SG discretization, HFS law, mobility
parameters, contact-row semantics and production fallback.  Its only change is
the geometry-defined transport-edge HFS drive

`sum(A_c * |grad(phin)_c|) / sum(A_c)`

in place of the existing magnitude of the area-weighted gradient vector.  No
target residual was used to choose an edge, branch, coefficient or support.

The candidate passes every aggregate and hotspot-migration protection check,
but it fails the mandatory per-node `<= 0.5x` gate at both endpoints.  At
`Vg=1.000000 V`, nodes 3721, 4091 and 3771 fail.  At `Vg=1.166667 V`, node
3721 remains just above the threshold at `0.507344x`.  Under the frozen stop
rule this is a negative T5 result: do not nominate a production C++ change,
do not enter T6, and do not run a reclose or 31-point curve.

## Pre-registration and exact reconstruction

`reference_staging/templates_ldmos_g3_wp3_t5_hfs_targeted_20260901/protocol.json`
was written before any candidate result was computed.  Its SHA-256 is
`6fa7f1ef32ffd0b82cf5caa5863cdeccd154b5eec928e981edbf2b88f6aead79`.
The protocol freezes one candidate, two endpoints, the T3 `0.5x` main gate,
all three migration-protection regions, and extreme-reconstruction tolerances
of relative L2 `<= 1e-10` and maximum physical error `<= 1e-18 A/um`.

Four same-state SG tables were generated: `edge_projection` and
`transport_cell_vector` at each endpoint.  Reconstructing their production
carrier residuals from edge incidence, including Ohmic carrier Dirichlet-row
replacement, gives:

| Vg (V) | Extreme | relative L2 | max abs (A/um) | Result |
|---:|---|---:|---:|---|
| 1.000000 | edge projection | 1.575726e-13 | 1.923875e-21 | pass |
| 1.000000 | transport-cell vector | 1.785680e-13 | 1.489874e-21 | pass |
| 1.166667 | edge projection | 2.006178e-13 | 5.264559e-21 | pass |
| 1.166667 | transport-cell vector | 2.117234e-13 | 6.527321e-21 | pass |

The independently reconstructed HFS formula agrees with the probe to maximum
relative error `4.439381e-16` and `4.413081e-16` at the two endpoints.  The
candidate reconstructed independently from each extreme agrees to
`1.398809e-14` and `1.415324e-14`.  Thus the negative gate result is not caused
by an algebraic-reconstruction mismatch.

## Frozen seven-node gate

| Vg (V) | seven-node L2 ratio | seven-node max ratio | nodes over 0.5x | Main gate |
|---:|---:|---:|---|---|
| 1.000000 | 0.215193 | 0.240745 | 3721: 1.186087; 4091: 0.796136; 3771: 1.999495 | fail |
| 1.166667 | 0.070646 | 0.053765 | 3721: 0.507344 | fail |

Complete per-node absolute ratios are:

| Vg (V) | 3721 | 3974 | 3973 | 4091 | 3727 | 4021 | 3771 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.000000 | 1.186087 | 0.069121 | 0.030735 | 0.796136 | 0.240745 | 0.168691 | 1.999495 |
| 1.166667 | 0.507344 | 0.094702 | 0.067529 | 0.249113 | 0.037245 | 0.081332 | 0.101972 |

The favorable seven-node aggregate ratios cannot override even one failing
point under the pre-registered conjunctive gate.

## Migration protection

| Vg (V) | candidate top-7 L2 / max | interface-band L2 / max | all-silicon L2 / max | top-7 overlap | Result |
|---:|---:|---:|---:|---:|---|
| 1.000000 | 0.999999462 / 0.999999848 | 0.951546 / 0.999999848 | 0.947413 / 0.999999848 | 7/7 | pass |
| 1.166667 | 0.999998676 / 0.999999626 | 0.842959 / 0.999999626 | 0.840547 / 0.999999626 | 5/7 | pass |

All candidate-own top-seven, interface-band and all-silicon L2/max metrics are
non-worsening.  At the upper endpoint, nodes 4676/4677 leave the top seven and
699/4569 enter it; the protected maxima still decrease slightly.  This passes
the migration guard but does not rescue the failed main gate.

## Implementation and artifacts

`scripts/audit_templates_ldmos_g3_wp3_t5_hfs_targeted.py` performs the
diagnostic-only reconstruction.  It verifies frozen input hashes, generates
the four SG probes, reproduces both registered extreme residuals, computes the
magnitude-average HFS drive directly from the sealed mesh and fixed state,
reconstructs the candidate from both extremes, and applies the closed gate.
It does not modify production C++ or defaults.

Primary outputs are under
`reference_staging/templates_ldmos_g3_wp3_t5_hfs_targeted_20260901/`:

- `protocol.json`: pre-result machine-readable protocol;
- `probes/<endpoint>/<extreme>/sg_edges.csv`: four SG edge tables, with sealed
  configs and stdout/stderr beside them;
- `analysis/extreme_reconstruction.csv`: exact-limit validation;
- `analysis/candidate_edges.csv` and `analysis/candidate_nodes.csv`: audit
  surfaces;
- `analysis/summary.json`: formal verdict and all gate metrics.

`tests/regression/test_audit_templates_ldmos_g3_wp3_t5_hfs_targeted.py`
covers the magnitude average, HFS reconstruction, carrier Dirichlet-row
semantics and the rule that aggregate improvement cannot hide a failing
point.  The T5 regression's four tests pass.  The final machine-readable
verdict is `stop_t5_h1_targeted_control_no_t6`.
