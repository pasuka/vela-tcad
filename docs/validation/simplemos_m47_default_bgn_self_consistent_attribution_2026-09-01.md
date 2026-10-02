# SimpleMOS M47 default-BGN self-consistent attribution

## Technical summary

The M46 peak is reproduced exactly at n23, Vd=0.05 V, Vg=0.05 V: Sentaurus gives 8.586250997030e-17 A/um and Vela gives 1.104643803536e-16 A/um, a 0.109418680924 dex (28.653%) Vela-over-Sentaurus difference.

The QF-referenced source barrier is -1.291129 mV lower in Vela. Its Maxwell-Boltzmann current proxy is only 0.021690 dex (19.8% of the observed log-current difference), leaving 0.087729 dex unattributed.

No potential, carrier-density, SRH, or quasi-Fermi metric satisfies the frozen target-localization rule. Their n23 discrepancies vary smoothly from Vg=0.00 to 0.10 V, while the current discrepancy has a sharp interior maximum at 0.05 V. The matched n19 control is also smooth. M47 therefore does not attribute the peak to a localized default-BGN-on self-consistent state error.

## The state differences form a smooth background, not the current peak

| Device | Vg (V) | Sentaurus Id (A/um) | Vela Id (A/um) | Current error (dex) | Barrier proxy (dex) | Integrated SRH ratio (dex) |
|---|---:|---:|---:|---:|---:|---:|
| n23 | 0.00 | 5.068684493e-17 | 5.273933704e-17 | 0.017239 | 0.020488 | -0.010026 |
| n23 | 0.05 | 8.586250997e-17 | 1.104643804e-16 | 0.109419 | 0.021690 | -0.009905 |
| n23 | 0.10 | 2.494951105e-16 | 2.730137802e-16 | 0.039123 | 0.022607 | -0.009772 |
| n19 | 0.00 | 1.729082366e-15 | 1.823920041e-15 | 0.023190 | 0.014508 | 0.004680 |
| n19 | 0.05 | 5.531688367e-15 | 5.862819933e-15 | 0.025249 | 0.014844 | 0.004656 |
| n19 | 0.10 | 1.821215226e-14 | 1.934788180e-14 | 0.026272 | 0.014969 | 0.004630 |

The barrier proxy increases monotonically across the three n23 gates (0.020488, 0.021690, and 0.022607 dex), so it cannot reproduce the 0.05 V current maximum. At the target, the barrier-node electron-density ratio is +0.021694 dex, but the mean over the source-barrier support is -0.013915 dex. This spatial sign change also argues against a single uniform carrier-density offset.

SRH is nearly identical spatially (Pearson r = 0.999817614) and its integrated magnitude differs by only -0.009905 dex. The electron and hole quasi-Fermi reference-plus-increment roundtrip is preserved to 6.939e-18 V.

The table is used instead of a trend chart because the contract contains only six discrete audit points and explicitly forbids interpolation.

## Contact flux shows endpoint redistribution, not a uniform scale error

At the target drain, the total-current error is 0.109419 dex and the electron component is 0.109456 dex; the hole component is negligible at this operating point. At the source, however, the total-current magnitude ratio is -0.102794 dex, in the opposite direction.

The signed source-plus-drain total is -3.251631860137e-17 A/um in Sentaurus and 1.703578187106e-17 A/um in Vela. This two-contact sum is an endpoint-partition diagnostic, not a KCL residual: gate/substrate terminals and volumetric generation are intentionally outside that sum. The opposite source/drain shifts show that the remaining 0.087729 dex cannot be described as a uniform source-to-drain flux multiplier.

## Scope, definitions, and method

M47 compares six exact self-consistent states: n23 and matched low-NWell n19 at Vd=0.05 V and Vg=0.00/0.05/0.10 V. n19 matches n23 in gate oxide time and LDD dose; NWell is the selected process-grid control. The current metric is signed log10(|Id,Vela|/|Id,Sentaurus|). The barrier proxy is -Delta(max(phin-psi))/(Vt ln(10)). Density and SRH ratios are signed log10 magnitude ratios.

Sentaurus T-2022.03-SP2 reused the M8 process TDRs and the complete original 0-to-2.5 V gate continuation, exporting only the three contracted gate states. Vela reused the frozen M8 Vd=0.05 V drain-ramp states and production configuration. Common Silicon node IDs and coordinates were compared directly without interpolation.

## Robustness and limits

All six Vela and Sentaurus currents reproduce the frozen M46 anchors within 1e-10 dex; all six states converge; coordinates are identical; required BGN, potential, density, SRH, and quasi-Fermi fields are present. The contract hash stayed unchanged. No HFS, SG, contact-current extraction, quasi-Fermi packing, BGN, SRH, mobility, solver, or production default was changed.

This is mechanism-consistent state attribution, not a causal intervention. A smooth field difference may contribute to the baseline cross-TCAD offset, but it does not explain why only the n23 0.05 V point rises to 0.109419 dex.

## Recommended next step

Keep the frozen defaults and closed topics intact. If a follow-on milestone is opened, isolate the 0.087729 dex residual through a fixed-state terminal-partition and local-continuity decomposition that includes the substrate terminal and volumetric source balance. Do not use it as authorization for HFS tuning or SG/contact extraction work.

## Further question

Does the n23 0.05 V residual arise from how the self-consistent bulk/substrate generation current is partitioned between source and drain, or from a transport response not captured by the scalar source-barrier proxy? M47 leaves that distinction unresolved.
