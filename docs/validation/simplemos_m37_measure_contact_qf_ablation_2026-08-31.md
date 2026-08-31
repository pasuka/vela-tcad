# SimpleMOS M37 measure, contact-SG, and QF-feedback ablation

## Scope

M37 separates three mechanisms at the n23 maximum-error point, Vd = 0.05 V
and Vg = 0.05 V.  The physics remains OldSlotboom + SRH(DopingDependence) +
PhuMob + Enormal + HighFieldSaturation(GradQuasiFermi), with
`transport_cell_vector` driving-field discretization.

The boundary-measure branch is a new self-consistent run.  The contact-SG and
quasi-Fermi branches reuse the accepted M22 and M26 frozen-operator/adjoint
evidence, so they retain exactly the production drain cut and current
functional.  No Sentaurus rerun or default-model change is involved.

## Independent ablations

| Mechanism | Intervention | Measured response | Result |
|---|---|---:|---|
| Boundary/contact measure | Signed measure only on external-boundary nodes and every contact one-ring | -0.000181 dex | 3.52% of the all-node signed-measure response |
| Contact SG flux | Fixed full state; HFS on/off on the native 7-edge drain cut | 0 A/um electron-cut change | All seven edge mobilities and electron fluxes unchanged |
| QF self-consistent feedback | M26 frozen direct term versus adjoint relaxation | 94.49% of finite HFS response | Dominant path; 5.50% first-order error |

For comparison, enabling signed AverageBox measure on all transport nodes
changes log current by `-0.0051444 dex`.  Restricting the replacement to the
external boundary and contact support changes it by only `-0.000180879 dex`.
Thus about 96.48% of the signed-measure response comes from nodes outside that
boundary/contact support.

The fixed-state contact-cut HFS replacement gives zero electron-current
change on every drain-cut edge.  The complete frozen terminal functional has
only a `-1.097e-21 A/um` response (`-4.31e-6 dex`), so this conclusion is not a
result of hiding a large direct contact contribution in a global sum.

The finite self-consistent HFS turn-on response is `-2.85476e-17 A/um`.  Its
direct frozen-operator contribution is only `-1.09768e-21 A/um`, or
`3.845e-5` of the finite response.  The adjoint relaxation term is
`-2.69751e-17 A/um`, explaining 94.49% within the measured 5.50% nonlinear
error envelope; 99.9997% of the signed adjoint contribution belongs to the
electron equation.

## Conclusion

The remaining SimpleMOS deep-off/HFS response is not controlled by an unknown
external-boundary measure or by direct HFS degradation on the drain-contact SG
edges.  The evidence instead supports an internal electron-continuity
perturbation that is transmitted to terminal current through self-consistent
quasi-Fermi relaxation.  This is the mechanism to test next in the LDMOS
case; the result does not justify changing a production default.
