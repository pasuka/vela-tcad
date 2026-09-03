# Genius NPN BJT transport and recombination acceptance

## Outcome

The 31-point M1 accepted-state chain passes convergence, three-terminal KCL,
terminal-current parity, and the potential/electron/hole spatial-state gates.
The expanded overall acceptance is **not yet passed** because the 3 V hole
current-density and Auger gates fail. SRH passes its magnitude, integral, and
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
- Vela authority: the 31/31 accepted-state manifest under `build-release/`
  (SHA-256 `0a1b36befc6be8bbebf728c861936e3aa82ac0fde3b46de38feaa43c1fee471b`).
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
- the Sentaurus-style nodal total current-density reconstruction.

All four exports contain 5611 common-mesh nodes. SRH+Auger closure and
drift+diffusion closure are exact at the precision of the exported values.

| VCE (V) | Vela total recombination integral (A/um) | SDevice comparison |
|---:|---:|---|
| 0 | 3.18029e-8 | available |
| 1 | 1.63754e-8 | available |
| 2 | 1.63657e-8 | missing SDevice state |
| 3 | 1.63574e-8 | available and asserted |

The existing SDevice run stored states at 0, 0.1, 1, and 3 V. The source deck
has been changed to store 0, 1, 2, and 3 V on the next trusted VM run. Current
SSH authentication rejects the documented private key, so no new 2 V oracle
was generated in this task.

## 3 V asserted comparison

| Quantity | Main observations | Gate |
|---|---|---:|
| Electron current density | P95 log error 0.27079 decade; normalized vector RMSE 0.56436; cosine 0.89655 | pass |
| Hole current density | P95 log error 1.48371 decades; normalized vector RMSE 0.81503; cosine 0.78774 | fail |
| SRH | integral ratio 0.77866; normalized L1 0.22140; peak offset 0.171875 um (diagnostic only) | pass |
| Auger | integral ratio 0.57542; normalized L1 0.42458; peak offset 1.59406 um | fail |

At 0 V, the P95 current-density magnitude errors are 0.06905 decade
(electron) and 0.08265 decade (hole). At 1 V they are 0.21464 and 0.50063
decade. The rapid increase of the hole-current discrepancy with collector bias
points to high-field/current-reconstruction or carrier-statistics spatial
alignment rather than a bias-independent unit conversion error.

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

## Remaining diagnosis

A node-level audit resolves the misleading 3 V Auger peak-location result. At
the SDevice argmax (node 1478), SDevice/Vela electron densities are
`5.3552900e19/5.3552910e19 cm^-3`, hole densities are
`1.2776783e14/1.2732199e14 cm^-3`, but Auger rates are
`1.0626467e23/5.0705278e22 cm^-3 s^-1`. At the Vela argmax (node 3549), the
same pattern remains: carrier densities agree within about 0.4%, while Auger
rates are `1.0427370e23/5.2135148e22 cm^-3 s^-1`. The SDevice peak values at
the two distant nodes differ by less than 2%, so a single-node argmax is
unstable and is retained only as a diagnostic, not an acceptance gate.

Using the configured electron Auger coefficient and the classical
`n*p-ni_eff^2` excess product predicts approximately the SDevice rate at these
high-electron-density nodes. Source inspection shows that Vela instead feeds
the Fermi-generalized SRH excess product into Auger. The evidence therefore
strongly supports an Auger excess-product/formulation mismatch as the dominant
amplitude error; it does not support a carrier-hotspot displacement or a simple
coefficient change. The next implementation task is an opt-in classical-Auger
A/B run with a matching analytic-Jacobian audit before changing any default.

The hole-current comparison should next audit the nodal reconstruction against
edge Scharfetter-Gummel fluxes in the base-collector depletion region. Vela's
drift/diffusion diagnostics are available but do not yet carry a verified
physical normalization, and the retained SDevice TDRs contain only total
carrier current density. A component-level cross-tool assertion therefore
requires both a verified Vela normalization and an SDevice decomposition or a
documented reconstruction.
