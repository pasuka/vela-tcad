# Templates/LDMOS G3 WP3 T4 physical-knockout matrix

Date: 2026-09-01
Frozen plan: `01f20ace88dc2b68a7bc3788534244dadf0d18f2` (R4)
Status: complete; all three branches are contract-not-equivalent

## Outcome

K1 HFS-off, K2 Boltzmann and K3 BGN-off each completed one independent
T-2022.03-SP2 two-endpoint solve on the sealed phase-01 grid.  All six TDRs
were imported independently and replayed only with the matching Vela physics.
Every endpoint fails at least one pre-registered cross-engine qualification
gate.  Therefore T4 provides **no authorized carrier/rule-out verdict** for
H1 or H2: the residual-reduction numbers below are diagnostic observations
only and must not be used for attribution or washout.

| Variant | Vg (V) | n/p P95 (dex) | n/p max (dex) | contact QF max (V) | neutrality max | port rel. | mesh/edges | residual reduction | R4 verdict |
|---|---:|---:|---:|---:|---:|---:|---|---:|---|
| K1 HFS-off | 1.000000 | 1.890710e-3 | 2.491747e-3 | 4.16e-16 | 1.01e-14 | 5.970707e-4 | exact | 5.73427x | contract-not-equivalent |
| K1 HFS-off | 1.166667 | 1.890711e-3 | 2.491748e-3 | 4.16e-16 | 1.01e-14 | 5.970850e-4 | exact | 17.83096x | contract-not-equivalent |
| K2 Boltzmann | 1.000000 | 4.817971e-1 | 1.319711 | 4.02e-16 | 1.64e-14 | 5.347107 | exact | 0.897530x | contract-not-equivalent |
| K2 Boltzmann | 1.166667 | 4.817971e-1 | 1.319711 | 4.02e-16 | 1.64e-14 | 5.347667 | exact | 1.012569x | contract-not-equivalent |
| K3 BGN-off | 1.000000 | 1.448334 | 2.902193 | 4.16e-16 | 1.38e-14 | 0.815040 | exact | 1.048315x | contract-not-equivalent |
| K3 BGN-off | 1.166667 | 1.448334 | 2.902193 | 4.16e-16 | 1.38e-14 | 0.815042 | exact | 1.010301x | contract-not-equivalent |

The R4 limits are P95 `<=1e-3 dex`, max `<=1e-2 dex`, contact QF
`<=1e-9 V`, neutrality `<=1e-6`, port relative error `<=1e-3`, and exact
mesh/transport-edge identity.  K1 misses only the n/p P95 gate.  K2 and K3
miss both n/p reconstruction and fixed-state port-current gates by large
margins.  Contact rows and topology pass in every case.

## VM execution and sealing

| Variant | Output name | host wall (s) | TDR count | deck SHA-256 | endpoint TDR SHA-256 (1.0 / 1.166667 V) |
|---|---|---:|---:|---|---|
| K1 | `g3_wp3_t4_k1_hfs_off_20260901_205717` | 262.94 | 4 | `9b727ff8...5c0d` | `7ab99e44...46aa` / `830eaac0...9d6` |
| K2 | `g3_wp3_t4_k2_boltzmann_20260901_210500` | 269.99 | 4 | `753acc4c...ef0` | `b39b1562...6c4e` / `00861a81...cce5` |
| K3 | `g3_wp3_t4_k3_bgn_off_20260901_211100` | 258.76 | 4 | `cc4f7326...8d3` | `977d96ed...0228` / `cb292621...4f34` |

All three remote runs returned zero, with no deck or license error.  The raw
archives contain the uploaded isolated deck and parameter file, a remote copy
of the already sealed `n1_fps.tdr`, both endpoint TDRs, PLTs, logs, timing and
exit code.  The grid SHA-256 is `07af6353c47c03a740583b873bac58e733a377a8bb9fa7e08baa3d30bc2c67f3`
in every run.

## Replay and gate implementation

`scripts/analyze_templates_ldmos_g3_wp3_knockout.py` performs the previously
missing TDR-to-R4 path:

1. verify the completed capture and prepared state-deck manifest hashes;
2. select exactly one TDR for each registered endpoint and import it with the
   frozen coordinate/doping policy;
3. build the complete Sentaurus `psi/phin/phip/n/p` state and run Vela
   `newton_carrier_term_probe` plus `sg_edge_flux_probe` under the matching
   K1/K2/K3 physics;
4. evaluate free-silicon n/p reconstruction, converged Ohmic boundary-state
   QF/neutrality, drain-cut port current, exact coordinate/edge signatures,
   fixed-seven and variant-top-seven residuals, and the frozen 40-edge
   denominator;
5. call the pre-registered qualifier and assert every produced summary against
   the closed R4 schema surface (including const/enum fields and artifact
   SHA-256 maps).

Ohmic neutrality uses the independently converged endpoint boundary `n/p`.
The SG probe's endpoint densities are bulk-kernel reconstructions and do not
apply Vela's Ohmic Dirichlet neutrality reconstruction; using them for the
contact-row gate would incorrectly mix gate (a) into gate (b).

The six formal summaries are:

- K1: `reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/g3_wp3_t4_k1_hfs_off_20260901_205717/t4_analysis_v4/{vg1p0,vg1p166667}/analysis/summary.json`
- K2: `reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/g3_wp3_t4_k2_boltzmann_20260901_210500/t4_analysis_v1/{vg1p0,vg1p166667}/analysis/summary.json`
- K3: `reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/g3_wp3_t4_k3_bgn_off_20260901_211100/t4_analysis_v1/{vg1p0,vg1p166667}/analysis/summary.json`

All six summaries pass the pre-registered
`vela.templates_ldmos.g3_wp3_knockout_summary.v1` schema contract.  The
analyzer regression plus qualifier/preparer regression totals 10 passing
tests.  IALMob, predictor, reclose, 31-point IdVg and production-default
changes remained disabled throughout.

## Decision

T4 is complete as an execution task, but inconclusive as a causal knockout
matrix because every branch failed the mandatory physical-equivalence
qualification.  In particular, K1's apparently large residual reductions
cannot promote H1, and K2/K3's near-baseline residual ratios cannot rule out
H2.  The next causal decision must use the independent T2 structure together
with a contract-equivalent targeted T5 control; no T4 branch is eligible for
direct reclose or 31-point promotion.
