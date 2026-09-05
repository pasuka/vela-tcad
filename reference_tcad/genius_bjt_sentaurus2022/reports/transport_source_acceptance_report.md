# Genius NPN BJT transport and recombination acceptance

## Outcome

The 31-point M1 accepted-state chain passes convergence, three-terminal KCL,
terminal-current parity, and the potential/electron/hole spatial-state gates.
The expanded overall acceptance is **not yet passed** because the 3 V hole
current-density gate fails. Both SRH and Auger now pass their integral and
normalized-shape gates.

The transport/source thresholds were introduced after the original
characterization. They are initial engineering regression gates rather than a
blind independent validation.

## Data provenance

- SDevice authority: the accepted T-2022.03-SP2 WP0-WP2 run recorded by
  `manifests/wp0_wp2_manifest.json` (SHA-256
  `9a5384c8e6cf492e07f983380002d5962191d658fa52b1aa8393f8b1ddfc6999`).
- Imported 3 V common-mesh field manifest: SHA-256
  `f728ffb14003d9d81840d4f2f99a07cf2e9a4c8ec8ebc857e6be9e20df115159`.
- Vela authority: the classical-Auger 31/31 accepted-state manifest under
  `build-release/` (SHA-256
  `7fc5286a3a0edc2e0449394dc129a9c4ac6e082557917051c5838a20672a43dd`).
- Comparison contract: schema 4; all source paths and hashes are repeated in
  the generated JSON summaries.

## Formal accepted-state export

`newton_solve_from_state` previously selected the basic VTK writer and omitted
the solver's existing transport/source diagnostics. It now uses the full
physics-aware writer. Each representative Vela state at VCE=0, 1, 2, and 3 V
exports:

- SRH, Auger, and total recombination rate;
- the corresponding signed electron and hole continuity source contribution;
- electron and hole drift, diffusion, and total vectors in Vela's legacy
  diagnostic scale, for algebraic closure only;
- electron, hole, and total current-density vectors recovered directly from
  the production SG line fluxes by a dual-face-weighted nodal least-squares
  projection. The SG line fluxes remain the conservative authority; these
  node vectors are comparison diagnostics.

All four exports contain 5611 common-mesh nodes. SRH+Auger closure and
drift+diffusion closure are exact at the precision of the exported values.

| VCE (V) | Vela total recombination integral (A/um) | SDevice comparison |
|---:|---:|---|
| 0 | 3.29434e-8 | available |
| 1 | 1.74925e-8 | available |
| 2 | 1.74825e-8 | missing SDevice state |
| 3 | 1.74740e-8 | available and asserted |

The existing SDevice run stored states at 0, 0.1, 1, and 3 V. The source deck
has been changed to store 0, 1, 2, and 3 V on the next trusted VM run. Current
SSH authentication rejects the documented private key, so no new 2 V oracle
was generated in this task.

## 3 V asserted comparison

| Quantity | Main observations | Gate |
|---|---|---:|
| Electron current density | P95 log error 0.01577 decade; normalized vector RMSE 0.03515; cosine 0.99939 | pass |
| Hole current density | P95 log error 0.82799 decade; normalized vector RMSE 0.06221; cosine 0.99810 | fail (P95 only) |
| SRH | integral ratio 0.77866; normalized L1 0.22140; peak offset 0.171875 um (diagnostic only) | pass |
| Auger | integral ratio 0.98781; normalized L1 0.01219; shape TV 0.00658; peak offset 0 | pass |

At 0 V, the P95 current-density magnitude errors are 0.01778 decade
(electron) and 0.01642 decade (hole). At 1 V they are 0.01445 and 0.24744
decade. At 3 V, the electron field and the hole-field normalized RMSE/cosine
gates pass; only the hole log-magnitude P95 gate remains. Its increase with
collector bias is localized to the low-current tail rather than a global
direction or scale mismatch.

## Nc/Nv and BGN-Fermi A/B result

The original Vela material used `Nc=2.8e19 cm^-3` and `Nv=1.04e19 cm^-3`.
Aligning these to the T-2022.03-SP2 values (`2.8567e19` and `3.1046e19 cm^-3`)
changed the 3 V base current from 6.33285e-8 to 5.98970e-8 A/um, versus the
SDevice value 6.01876e-8 A/um. The SRH and Auger integral ratios improved from
0.76815 to 0.77866 and from 0.56002 to 0.57542. This alignment was adopted and
the full 31-point chain was rerun.

With aligned Nc/Nv, disabling the OldSlotboom Fermi-statistics correction
reduced the base current to 4.06488e-8 A/um. Although some source-distribution
metrics improved, terminal parity degraded materially; the correction therefore
remains enabled.

## Classical Auger A/B and Jacobian audit

A node-level audit resolves the misleading 3 V Auger peak-location result. At
the SDevice argmax (node 1478), SDevice/Vela electron densities are
`5.3552900e19/5.3552910e19 cm^-3`, hole densities are
`1.2776783e14/1.2732199e14 cm^-3`, but Auger rates are
`1.0626467e23/5.0705278e22 cm^-3 s^-1`. At the Vela argmax (node 3549), the
same pattern remains: carrier densities agree within about 0.4%, while Auger
rates are `1.0427370e23/5.2135148e22 cm^-3 s^-1`. The SDevice peak values at
the two distant nodes differ by less than 2%, so a single-node argmax is
unstable and is retained only as a diagnostic, not an acceptance gate.

The implementation now exposes `auger_excess_product`. Its compatibility
default is `generalized_fermi`; this BJT M1 fixture opts into `classical_np`.
At 3 V this changes the Auger integral ratio from 0.57542 to 0.98781, normalized
L1 from 0.42458 to 0.01219, shape TV from 0.09196 to 0.00658, and aligns the
reported peak to node 1478. The full 31-point chain remains accepted; 3 V
collector/base currents change by less than 2.4e-15/4.3e-17 A/um relative to
the previous accepted state.

The analytic Jacobian uses `p*dn+n*dp` for the classical excess product. A
source-only centered-difference audit on the real 3 V state at nodes 1478 and
3549 gives an overall relative difference of `5.07e-9`; the psi, phin, and
phip column-block differences are `2.32e-9`, `5.92e-9`, and `3.52e-9`.
The small-mesh Fermi/Auger regression test independently passes.

The unweighted Auger node-log P95 remains 0.903 decade because nodes down to
1e-6 of the reference peak are counted equally; at the 1% peak mask it is
0.215 decade. Since the absolute integral differs by only 1.22% and normalized
shape error is 0.0122, node-log P95 and single-node argmax are retained as
diagnostics rather than source acceptance gates.

## Hole-current edge audit

The base-collector window audit selects 516 vertically oriented edges at
`x=2.5..4.5 um`, `y=0.55..0.85 um`; 43 cross the net-doping sign change. On
the 86 edges above 1e-3 of the SDevice projected-current peak:

- Vela nodal reconstruction versus SDevice nodal projection has normalized
  RMSE 0.06276, cosine 0.99925, and fitted scale 1.05194.
- Production SG edge flux versus SDevice nodal projection has normalized RMSE
  0.71892, cosine 0.99675, and requires a fitted scale of 3.5301.
- Production SG flux versus Vela nodal reconstruction has normalized RMSE
  0.70322, cosine 0.99889, and fitted scale 3.3605.

This evidence motivated the implemented recovery. Each production SG line
flux is divided by its dual-face length to obtain an edge-normal current
density, then the incident edge values are fitted at each node using the same
dual-face lengths as least-squares weights. An independent Python replay and
the C++ VTK output agree to `1.88e-13 A/cm^2` maximum absolute difference.

At 3 V, the new recovery reduces electron normalized RMSE from 0.56436 to
0.03515 and hole normalized RMSE from 0.81503 to 0.06221; cosine similarities
rise to 0.99939 and 0.99810. The pre-registered hole-current mask and 0.5-decade
P95 limit were deliberately retained. Hole P95 is still 0.82799 decade at the
`1e-6` peak-relative mask, but falls to 0.35064 at `1e-4` and 0.04269 at `1e-2`.
The remaining failure is therefore a low-current-tail magnitude discrepancy,
not a principal-current-vector recovery failure. SDevice exposes only nodal
current vectors here, so its internal directed-edge flux remains unavailable
as a strict edge oracle.

## Low-current tail localization

The unchanged `1e-6` mask selects 2121 nodes; 180 nodes (8.49%) exceed
0.5 decade. All 180 have lower Vela magnitude and lie in the p-type base; 179
are interior nodes, and 144 lie at `y=0.50..0.75 um`. By SDevice current
direction, 141 are `+x`, 27 are `-y`, and 12 are `+y`. No n-type emitter-side,
contact, contact-one-ring, or emitter-base-junction node exceeds 0.5 decade.

The discrepancies are explicitly confined to weak currents. All 1700 nodes at
or above `1e-2` of the SDevice peak stay below 0.5 decade and contain 99.991%
of the reference-current squared magnitude. In the `[1e-5,1e-4)` and
`[1e-6,1e-5)` peak-relative bins, 53/63 and 50/52 nodes exceed 0.5 decade,
respectively.

Three checks distinguish the remaining hypotheses:

- Five nodal recovery weights change global P95 by only 0.0069 decade, so the
  choice among uniform/primal/dual weighting does not explain the tail.
- On the 180 tail nodes, SDevice/Vela hole density, mobility, and common-mesh
  quasi-Fermi-gradient P95 differences are only 0.0322, 0.0043, and 0.0070
  decade. These local model factors are too close to explain the median
  1.241-decade current ratio directly.
- SDevice current on these nodes is a median 22.2 times its own
  `q*mu*p*|grad(phi_p)|` nodal proxy, while Vela is 1.000 times its proxy. A
  SDevice nodal-vector projection/recovery round trip has only 0.0461-decade
  global P95 and does not make the Vela comparison closer.

The strongest current evidence therefore favors an SDevice-specific
weak-current element-to-node construction or extrapolation effect, rather than a
Vela local-state, mobility, or SG recovery-weight defect. This remains a
best-supported working hypothesis rather than a strict proof because the stored
SDevice TDR does not expose its directed-edge transport flux. The formal gate
is intentionally unchanged and remains failed. Reproduce with
`scripts/diagnose_genius_bjt_hole_current_tail.py`; full strata and source
hashes are recorded in `reports/hole_current_tail_diagnosis.json`.
