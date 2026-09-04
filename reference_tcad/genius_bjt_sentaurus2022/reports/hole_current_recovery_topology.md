# Genius NPN BJT weak-hole-current recovery topology audit

## Technical summary

The coarse-grid 0.828-decade P95 failure is produced after the conservative SG edge fluxes are mapped to a vertex vector. Replaying the production dual-face least-squares formula from the saved edges reproduces the exported Vela field to 1.847e-13 A/cm2. Changing only the recovery order to per-cell fitting followed by area-weighted cell-to-node projection reduces the original 180-node tail P95 from 1.869552 to 0.121906 decade. The SDevice node source used here matches the formal comparison source with normalized vector RMSE 4.797e-11.

## Recovery-order comparison

| Population | Nodes | Direct nodal P95 (dec) | Cell-first P95 (dec) | Direct failures | Cell-first failures |
|---|---:|---:|---:|---:|---:|
| dominant_current | 1700 | 0.042748 | 0.023536 | 0 | 0 |
| weak_main_current | 421 | 1.663415 | 0.114394 | 180 | 0 |
| original_tail | 180 | 1.869552 | 0.121906 | 180 | 0 |

The dominant-current region is insensitive to recovery order. The improvement is concentrated in the weak-current base. The cell-first path reconstructs each triangle from all three of its edges, including the edge opposite the target vertex, before projecting the cell vectors to that vertex; the direct-node path uses only edges incident on the vertex.

## Topology evidence

| Population | Direct-fit residual median | Cell cancellation median | Cell-vector dispersion median | Fit condition median |
|---|---:|---:|---:|---:|
| dominant_current | 0.039926 | 1.00098 | 0.108259 | 1.5 |
| weak_main_current | 0.151124 | 1.00274 | 1.49481 | 1.5 |
| original_tail | 0.107865 | 1.00009 | 1.84691 | 2 |

The condition-number comparison tests mesh-direction degeneracy; the residual and cell-dispersion comparisons test whether one smooth vector can represent the whole vertex patch. Tail cell-vector dispersion is much larger than in the dominant region, but the weak within-tail correlations mean it classifies the affected region better than it predicts node-by-node error severity. The near-unity cancellation factor and modest condition numbers reject simple vector cancellation and ill-conditioning as primary explanations.

## Correlations on the original tail

- direct_error_vs_log10_cell_cancellation: 0.102238
- direct_error_vs_log10_cell_vector_dispersion: 0.082660
- direct_error_vs_log10_direct_fit_residual: 0.216787
- direct_error_vs_log10_fit_condition_number: -0.399711

## Assessment

- Verified: the formal Vela node field is an exact replay of the saved SG edge projections under the production direct-node formula.
- Verified: using the same edge fluxes, cell-first recovery removes most of the weak-current discrepancy without changing the state, mobility, or flux operator.
- Strongly supported: the remaining coarse-grid P95 failure is a representation-support mismatch between Vela's direct vertex fit and SDevice's element-oriented current field, amplified where neighbouring cell currents vary in direction.
- Not proven: the exact proprietary SDevice element-to-vertex weighting. Terminal currents and conservative section fluxes remain the physical acceptance oracle.
