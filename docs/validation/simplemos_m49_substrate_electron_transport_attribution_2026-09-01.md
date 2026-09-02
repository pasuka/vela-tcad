# SimpleMOS M49 substrate electron-transport attribution

## Technical summary

M49 reuses the six exact M47 default-BGN-on states and the frozen M48 substrate-electron terminal anchors. The outcome is `boundary_observable_limited`. The largest exact bridge term is `terminal_boundary_residual` at 99.99% of the target terminal difference.

The direct boundary-field reconstruction leaves 99.99% of the target terminal difference in the terminal-to-boundary residual. Therefore the analysis reports where the fixed-state discrepancy is visible without converting a field-correlation result into a code-defect claim.

Only 0.0058% of the target terminal difference appears in the cross-solver boundary integral. Across the three n23 gate points that boundary-field difference varies by just 0.067%, while the target terminal difference is 278.3x the Vg=0 value and 4.6x the Vg=0.1 value. The exported boundary field therefore does not carry the target localization.

The layer-0 SRH magnitude ratio is floor-sensitive because the p95 absolute rates are only 6.978410e-05 cm^-3 s^-1 in Sentaurus and 0.000000e+00 cm^-3 s^-1 in Vela. It is not used as a mechanism discriminator; M48 already established that the integrated SRH difference is nonmaterial.

## Exact target bridge

| Term | Contribution (A/um) | Absolute fraction |
|---|---:|---:|
| density | 3.628404986888e-24 | 0.000 |
| mobility | -1.030277832265e-20 | 0.000 |
| qf_gradient | -9.230044756815e-21 | 0.000 |
| constitutive_residual | 1.679547713847e-20 | 0.000 |
| terminal_boundary_residual | -4.694430122779e-17 | 1.000 |

The five terms sum to -4.694703494533e-17 A/um, matching the M48 terminal difference -4.694703494533e-17 A/um with relative error 0.000e+00.

## Substrate graph-distance rows

| Layer | Nodes | QF-gradient p95 delta (V/cm) | n p95 (dex) | mobility p95 (dex) | Jn p95 (dex) | psi p95 (mV) | SRH p95 (dex) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 29 | 6.116177e+03 | 0.000006 | 0.107120 | 1.026320 | 0.000119 | 11.964746 |
| 1 | 29 | 6.113435e+03 | 0.000277 | 0.105019 | 0.979422 | 0.000119 | 0.001464 |
| 2 | 29 | 6.122187e+03 | 0.000111 | 0.102904 | 0.899777 | 0.000119 | 0.000977 |
| 3 | 29 | 6.136234e+03 | 0.000127 | 0.101147 | 0.844029 | 0.000121 | 0.001016 |
| 4 | 29 | 6.160564e+03 | 0.000422 | 0.098634 | 0.711096 | 0.000130 | 0.001301 |

## Exact-state controls

| Device | Vg (V) | Terminal difference (A/um) | Dominant term | Dominant absolute fraction | Boundary residual fraction |
|---|---:|---:|---|---:|---:|
| n23 | 0.00 | -1.686669943551e-19 | terminal_boundary_residual | 0.984 | 0.984 |
| n23 | 0.05 | -4.694703494533e-17 | terminal_boundary_residual | 1.000 | 1.000 |
| n23 | 0.10 | -1.018341161792e-17 | terminal_boundary_residual | 1.000 | 1.000 |
| n19 | 0.00 | 3.764463040765e-20 | constitutive_residual | 0.871 | 0.323 |
| n19 | 0.05 | -1.698460041140e-17 | terminal_boundary_residual | 1.001 | 1.001 |
| n19 | 0.10 | -4.638204429937e-18 | terminal_boundary_residual | 1.005 | 1.005 |

## Scope and methodology

The substrate boundary is the set of Silicon boundary edges whose endpoints are both in the M47 substrate-contact node set. The outward normal is determined from the adjacent Silicon triangle. Exported nodal current density is integrated with endpoint-average edge quadrature.

The reconstructed conventional electron current uses Jn = q n mu grad(Phi_n). An exact three-factor Shapley bridge separates density, mobility, and normal quasi-Fermi-gradient. Two explicit residuals preserve the exported-field constitutive mismatch and the difference between boundary-field quadrature and the frozen M48 terminal observable.

## Robustness and claim boundary

All six cross-solver states use identical node IDs and coordinates; no interpolation is used. The factor bridge is algebraically exact. A large terminal-to-boundary residual triggers the contract's boundary-observable-limited outcome rather than being silently assigned to a transport factor.

M49 does not rerun either solver, change BGN/SRH/mobility/contact/HFS/default settings, or reopen SG, contact extraction, or quasi-Fermi packing. Fixed-state localization is descriptive and does not by itself identify a production-code defect.

## Recommended next step

If the boundary-observable gate passes, the dominant localized state factor is the appropriate candidate for a separately contracted intervention. If it fails, the next admissible task is a diagnostic-output replay that exports native substrate-face flux on the unchanged states; it is not a renewed contact-current extraction investigation.
