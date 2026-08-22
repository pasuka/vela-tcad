# Sentaurus 2022 workflow validation report

Validation date: 2026-08-22

Machine-readable closeout: `docs/validation/sentaurus2022_remaining_tasks_20260822.json`.

## Build and regression status

- Windows MSYS2 UCRT64 Release targets were rebuilt and used for all new
  full-mesh simulations.
- Python regression coverage includes frontend provenance, ABA/BAB geometry,
  contact edge splitting, ordered doping replacement, directed Gaussian
  profiles, refinement controls, execution-stage ordering, IV/BV mode mapping,
  contract loading, and machine-readable failure reports.

## Contract executions

| Contract case | Runner | Result | Evidence |
| --- | --- | --- | --- |
| `pn2d_sentaurus2018_coarse7x3_iv` | Debug (earlier baseline) | accepted | 21 nodes, 24 cells, 16 points, last bias 0.3 V, maximum terminal-current imbalance 0.001216 |
| `pn2d_sentaurus2018_coarse7x3_bv` | Debug (earlier baseline) | accepted | reached the contracted -0.2 V stop with `bv_reverse` output mode |
| `pn2d_sentaurus2018_iv` | Release | accepted | 1767 nodes, 3360 cells, 18 accepted points, last bias 0.3 V; maximum absolute terminal mismatch 1.637e-18 A/um and maximum mixed-tolerance ratio 0.5511 |
| `pn2d_sentaurus2018_bv` | Release | failed acceptance after reaching target | 1767 nodes, 3360 cells, 59 accepted solver points, last bias -20 V with `bv_reverse`; maximum absolute terminal mismatch 6.474e-4 A/um and mixed-tolerance ratio 20.0 |

The full PN2D IV low-current false rejection was caused by applying a pure
relative terminal-current ratio near zero current. Acceptance now retains the
raw relative diagnostic but applies a checked-in mixed absolute/relative gate.
The absolute allowance is 2e-18 A/um; the measured maximum is 1.637e-18 A/um.
Gummel carrier convergence is componentwise for this IV import, preventing the
global majority-carrier L2 norm from hiding a local minority-carrier update.

The full PN2D BV contract is now checked in and the solver reaches its -20 V
target. It is deliberately not marked `accepted`: at high reverse bias the
anode current is about 6.1e-4 A/um while the cathode current remains near zero,
so the terminal mismatch is physical/numerical rather than a near-zero
roundoff artifact. Failure reports preserve the last converged bias and remove
stale success manifests.

## Sentaurus 2022 official SDE oracle

The configured licensed VM ran Sentaurus Structure Editor T-2022.03-SP2. The
frozen `reference_tcad/pn2d_sentaurus2022_oracle` fixture contains a
variable-expanded `.log.eval`, official `_bnd.tdr`, `_msh.cmd`, `_msh.tdr`,
stdout/stderr, version output, expected Device IR and SHA-256 manifest. The
official mesh has 1943 points, 3680 elements and 3 regions. The VM evidence also
records T-2022.03's deprecation warning for the source deck's legacy
`sdegeo:define-2d-contact` calls.
