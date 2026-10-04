# SimpleMOS M48 terminal partition and continuity closure

## Technical summary

M48 exactly decomposes the M47 target drain-current gap of 2.460187038329e-17 A/um without running either solver again. The largest bookkeeping term is substrate rebalancing: 4.654871563240e-17 A/um, or 189.2% of the signed gap. It is offset by source rebalancing (-101.4%); the cross-solver KCL-residual difference contributes 12.2%.

The substrate term is target-localized under the frozen adjacent-gate and n19 rule: True. Both target states pass Id/|KCL| >= 10 (Sentaurus 28.589; Vela 3.660e+07). Therefore numerical nonclosure is measurable but does not dominate the gap.

The integrated SRH difference is only 1.57% of the drain gap and fails the materiality rule. M48 identifies a target-localized, electron-channel substrate/source redistribution in the terminal ledger; this is an exact continuity bookkeeping result, not proof of the physical operator that created the state.
The M47 0.087729 dex barrier-proxy residual is logarithmic and is not an additive term in this A/um identity. M48 localizes its terminal manifestation; it does not assign the whole 0.087729 dex to the substrate term.

## The target drain gap is balanced mainly by substrate and source redistribution

| Identity term | Contribution (A/um) | Signed fraction of drain gap |
|---|---:|---:|
| Substrate rebalancing | 4.654871563240e-17 | 1.892 |
| Source rebalancing | -2.495023008915e-17 | -1.014 |
| KCL residual difference | 3.003384840028e-18 | 0.122 |
| Gate rebalancing | -0.000000000000e+00 | -0.000 |

These four terms sum to the drain-current gap to the frozen 1e-12 relative identity tolerance. Fractions may exceed 100% because compensating signed terms are retained.

## The substrate contribution is localized at the M47 peak

| Device | Vg (V) | Drain gap (A/um) | Substrate/gap | Source/gap | KCL/gap | |Delta SRH|/|gap| |
|---|---:|---:|---:|---:|---:|---:|
| n23 | 0.00 | 2.052492e-18 | -0.113 | 1.122 | -0.008 | 0.189 |
| n23 | 0.05 | 2.460187e-17 | 1.892 | -1.014 | 0.122 | 0.016 |
| n23 | 0.10 | 2.351867e-17 | 0.416 | 0.542 | 0.042 | 0.016 |
| n19 | 0.00 | 9.483768e-17 | -0.009 | 1.009 | 0.000 | 0.002 |
| n19 | 0.05 | 3.311316e-16 | 0.049 | 0.951 | 0.000 | 0.001 |
| n19 | 0.10 | 1.135730e-15 | 0.003 | 0.997 | -0.000 | 0.000 |

The exact six-point table is used instead of an interpolated trend. The target substrate absolute gap fraction exceeds both adjacent n23 gates and matched n19 at Vg=0.05 V. The source term provides the compensating endpoint response. The KCL contribution is also target-localized by the same rule, but at 12.2% of the gap it remains a minority term.

## Carrier continuity separates substrate electron partition from SRH

At the target, Sentaurus electron continuity has a residual of -2.982332909964e-18 A/um and hole continuity has -2.105494846883e-20 A/um. Vela gives 9.350292078939e-21 and -9.353310484047e-21 A/um, respectively.

The substrate electron-current difference is -4.694703494533e-17 A/um, whereas the substrate hole-current difference is 3.983193129257e-19 A/um. The dominant partition therefore sits in the electron channel. The much smaller integrated SRH difference does not account for that redistribution.

## Scope, definitions, and method

M48 reuses the six exact M47 states for n23 and n19 at Vd=0.05 V and Vg=0.00/0.05/0.10 V. Sentaurus conventional total current is eCurrent+hCurrent. Vela stores the hole component with particle-current orientation, so its conventional hole current is the negative of the recorded hole column and total current is electron minus hole.

For each solver and state, electron closure is the sum of four conventional electron terminal currents plus integrated signed SRH current. Hole closure is the four-terminal conventional hole sum minus signed SRH current. The total identity is Delta Id = Delta KCL - Delta Isource - Delta Isubstrate - Delta Igate, with every Delta defined as Vela minus Sentaurus.

## Robustness, limitations, and claim boundary

All 144 terminal-component rows and 12 carrier-continuity rows are finite; all component and terminal identities pass; all drain currents reproduce M47 within 1e-10 dex; the contract hash stayed frozen. No solver, physics, extraction, or historical artifact changed.

A terminal identity determines where signed current is balanced, not why the nonlinear self-consistent solution chose that partition. The Sentaurus target KCL residual explains a minority of the cross-solver gap and must remain visible, but both solvers pass the contract's numerical-resolution gate.

Sentaurus carrier closure uses its exported nodal SRH field with the common mesh-volume quadrature, while Vela also exposes a native cell-integrated SRH balance. Thus the carrier-row closures are field-derived diagnostics, not a claim that Sentaurus's internal Newton residual has been observed. The four-terminal KCL identity does not depend on that SRH quadrature.

## Recommended next step

The terminal manifestation of the residual is now localized to the self-consistent electron-current partition between substrate, source, and drain. A further milestone, if desired, should compare the substrate-side electron quasi-Fermi/continuity rows on the same frozen states. It must not reopen SG extraction, contact extraction, quasi-Fermi packing, or global HFS tuning.

## Further question

Which substrate-side electron-continuity rows first create the target-localized partition difference, and is their response tied to BGN-dependent equilibrium state or to another fixed production operator? M48 does not make that causal assignment.
