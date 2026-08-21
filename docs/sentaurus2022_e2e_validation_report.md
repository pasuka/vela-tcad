# Sentaurus 2022 workflow validation report

Validation date: 2026-08-21

## Build and regression status

- Windows MSYS2 UCRT64 Debug build completed successfully.
- Full Debug CTest suite: 663/663 passed.
- An independent UCRT64 Release `vela_example_runner` build completed
  successfully and was used for the full-mesh simulation below.
- Python regression coverage includes frontend provenance, ABA/BAB geometry,
  contact edge splitting, ordered doping replacement, directed Gaussian
  profiles, refinement controls, execution-stage ordering, IV/BV mode mapping,
  contract loading, and machine-readable failure reports.

## Contract executions

| Contract case | Runner | Result | Evidence |
| --- | --- | --- | --- |
| `pn2d_sentaurus2018_coarse7x3_iv` | Debug | accepted | 21 nodes, 24 cells, 16 points, last bias 0.3 V, maximum terminal-current imbalance 0.001216 |
| `pn2d_sentaurus2018_coarse7x3_bv` | Debug | accepted | reached the contracted -0.2 V stop with `bv_reverse` output mode |
| `pn2d_sentaurus2018_iv` | Release | failed acceptance | 1767 nodes, 3360 cells, 16 solved points, last bias 0.3 V; the only failed gate was terminal-current imbalance (0.4373 versus 0.02 maximum) in the low-current range |

The full PN2D IV run is therefore a completed numerical execution, but it is
not reported as `accepted`. Tightening the Gummel iteration tolerance and
raising its iteration budget did not improve the low-current imbalance;
Newton was rejected and the solver used the validated Gummel fallback. Those
ineffective tuning changes were reverted. The checked-in acceptance threshold
was not relaxed.

The full PN2D BV case has no checked-in contract case and was not claimed as
complete. The licensed Sentaurus 2022 executable was not available on this
development machine, so official `.log.eval`, TDR, and mesh oracle generation
remains a VM validation step described in
`docs/sentaurus2022_vm_fixture_generation.md`.
