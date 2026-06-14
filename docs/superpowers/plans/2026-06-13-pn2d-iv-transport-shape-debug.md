# PN2D IV Transport Shape Debug Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Localize the remaining pn2d Sentaurus2018 IV current trough by comparing transport-driving physical quantities instead of tuning global model parameters.

**Architecture:** The electrostatic state is close enough that the next debug layer should compare derived transport quantities: current-density fields, quasi-Fermi gradients, contact-edge current concentration, and Scharfetter-Gummel flux terms. The plan first builds read-only diagnostics from existing Sentaurus exports and Vela fixed-bias VTK/CSV outputs, then uses one minimal perturbation at a time only if the diagnostics identify a specific failing layer.

**Tech Stack:** C++20, CMake/Ninja, Python regression diagnostics, CSV/JSON reports, Pillow plots, MSYS2 UCRT64 on Windows.

---

## Current Evidence and Priority

- `0 V` state quantities are largely aligned:
  - `NetDoping`: mean abs diff `15.7 cm^-3`, max `16 cm^-3`.
  - `Potential`: mean abs diff `8.24e-3 V`, p95 `9.85e-3 V`.
  - `ElectricField`: mean abs diff `434 V/cm`, p95 `1.40e3 V/cm`, max `2.93e3 V/cm`.
  - `eDensity/hDensity`: mean abs diff about `1.99e14 cm^-3`, p95 about `1.72e15 cm^-3`.
  - `electron_qf/hole_qf`: near numerical agreement, p95 below `4e-10 V`.
- Fixed same-bias `1 V` state quantities are also close enough to rule out gross electrostatic mismatch:
  - `Potential`: mean abs diff `2.57e-3 V`, max `9.85e-3 V`.
  - `electron_qf`: mean abs diff `7.61e-3 V`, p95 `1.31e-2 V`, max `1.53e-2 V`.
  - `hole_qf`: mean abs diff `7.69e-3 V`, p95 `1.31e-2 V`, max `1.51e-2 V`.
  - `ElectricField`: mean abs diff `104 V/cm`, p95 `333 V/cm`, max `802 V/cm`.
  - `eDensity/hDensity`: mean abs diff about `1.42e16 cm^-3`, p95 about `1.86e16 cm^-3`, around `8-11%` p95 relative error.
- The remaining large error is IV curve shape:
  - `0.25 V` Vela/Sentaurus ratio `0.998`.
  - `0.30 V` ratio `0.642`.
  - `0.80 V` ratio `0.557`.
  - `1.00 V` ratio `0.827`.
  - Bias bucket `0.30..0.80 V` mean ratio `0.442`, so this is the primary residual mismatch.
- Mobility and recombination/BGN scans did not identify a promotable single parameter:
  - `no_mobility` fixes the trough but overdrives `1 V` to `1.491x`.
  - SRH lifetime/BGN candidates either leave `0.8 V` unchanged or destroy the good `0.25 V` match.

## Execution Results

- Added `scripts/analyze_pn2d_iv_transport_shape.py`.
- Generated read-only reports under `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/`.
- Sentaurus `sim_fields/iv/fields` currently contains one region0 field state, treated as `source_bias_V=1.0`; it does not provide a multi-bias field sequence for `0.25/0.30/0.50/0.80 V`.
- Vela multi-bias proxy reports show:
  - `ElectronQuasiFermi` and `HoleQuasiFermi` spans track the applied bias smoothly.
  - Mean QF gradients increase monotonically from `0.25 V` to `1.0 V`.
  - Mean electron/hole mobility increases smoothly from `0.25 V` to `1.0 V`; no mobility trough matching the terminal-current ratio trough was observed.
  - Max electric field decreases smoothly with forward bias as the junction barrier collapses.
- Contact edge concentration:
  - Cathode `top3_abs_fraction` is stable near `0.1875` at `0.25/0.30/0.50/0.80/1.00 V`.
  - Contact edge sums match terminal balance at the audited points; no terminal-balance gate failures were found.
- Added opt-in C++ transport diagnostics columns:
  - `mean_electron_mobility_m2_V_s`
  - `mean_hole_mobility_m2_V_s`
  - `min_electron_mobility_m2_V_s`
  - `min_hole_mobility_m2_V_s`
  - `max_electric_field_V_per_cm`
  - `mean_electron_qf_gradient_V_per_cm`
  - `mean_hole_qf_gradient_V_per_cm`
- Decision: the next root-cause branch is `sg_flux_or_mobility_einstein_relation`. The current evidence continues to rule out contact-current extraction and nonlinear-state QF span as primary explanations.
- Follow-up execution of `sg_flux_or_mobility_einstein_relation`:
  - `ContactCurrent` and `CoupledDDAssembler` both use the variable-`ni` quasi-Fermi SG helpers; the coefficient chain is consistent (`mu * Vt / edgeLength`, with `edge.couple` applied once for contact integration).
  - Density-form SG and variable-`ni` quasi-Fermi SG agree on the audited contact edges at key biases: density/QF flux ratio is `1.000023` at `0.25 V`, `1.0000035` at `0.30 V`, and effectively `1.0` from `0.50 V` through `1.0 V`.
  - Vela contact current density at `1.0 V` is `20812.44 A/cm^2`; Sentaurus `TotalCurrentDensity` mean is `25175.98 A/cm^2`; ratio `0.826678` matches the terminal IV ratio, so width/current-density projection is not the residual error.
  - At `1.0 V`, Vela component current densities are low in both carriers: electron `14220.76 A/cm^2` versus Sentaurus `17494.42 A/cm^2` (`0.8129x`), hole `6591.67 A/cm^2` versus `7681.58 A/cm^2` (`0.8581x`).
  - Added generated reports:
    - `effective_mobility_proxy_1v.csv/json`
    - `contact_edge_effective_mobility_proxy_1v.csv/json`
    - `contact_edge_effective_mobility_proxy_1v_edges.csv`
  - On the majority-carrier extraction edges at `1.0 V`, Sentaurus has about `16-17%` larger `q * carrier_density * grad(QF)` than Vela:
    - Cathode electron mean `2.075e9` versus `1.786e9` in SI driver units.
    - Anode hole mean `2.076e9` versus `1.774e9`.
  - Majority-carrier effective mobility inferred from Sentaurus fields is close to or modestly above Vela's contact-current effective mobility:
    - Cathode electron Sentaurus median `0.07234 m^2/V/s`, Vela current-derived `0.06831 m^2/V/s`.
    - Anode hole Sentaurus median `0.03175 m^2/V/s`, Vela current-derived `0.03199 m^2/V/s`.
  - Current evidence points to a combined contact-edge transport-driver deficit in Vela: slightly lower contact-edge QF gradient/carrier driver plus a small mobility mismatch, not a sign convention, SG algebra, or current extraction bug.
- Next debug execution of the contact-edge transport-driver branch:
  - Added reproducible contact-edge proxy reports to `scripts/analyze_pn2d_iv_transport_shape.py`:
    - `contact_edge_transport_proxy_1v.csv/json`
    - `contact_edge_transport_proxy_1v_edges.csv`
    - `contact_edge_transport_proxy_compare_1v.csv/json`
  - Added regression coverage proving the script emits the contact-edge transport proxy comparison.
  - At `1.0 V`, the majority-carrier density on the extraction edges is almost aligned while QF edge drop is low in Vela:
    - Cathode electron `qf_drop_V`: Sentaurus mean `4.39326e-3 V`, Vela `3.81483e-3 V`, ratio `0.86834`.
    - Cathode electron density: ratio `0.99125`.
    - Cathode electron `q*n*grad(QF)`: ratio `0.86074`.
    - Anode hole `qf_drop_V`: Sentaurus mean `4.39500e-3 V`, Vela `3.78700e-3 V`, ratio `0.86166`.
    - Anode hole density: ratio `0.99181`.
    - Anode hole `q*p*grad(QF)`: ratio `0.85461`.
  - Node-level split shows the contact QF boundary values are aligned; the entire QF-drop deficit comes from the first interior nodes:
    - Cathode electron contact QF delta Vela-Sentaurus: approximately `0 V`; interior delta `-0.57843 mV`.
    - Anode hole contact QF delta Vela-Sentaurus: `0 V`; interior delta `+0.60800 mV`.
    - Majority carrier density at those first interior nodes is `1.5-1.6%` lower in Vela.
  - Local Sentaurus IV TDR evidence:
    - `pn2d_iv_des.tdr` inventory contains one final field state, with contact voltages `Cathode=0 V`, `Anode=1 V`.
    - `pn2d_iv_sdevice.cmd` writes a single `Plot = "pn2d_iv_des.tdr"` and does not use collected/multi-bias plot output.
    - Therefore field-level Sentaurus/Vela contact-edge transport comparison is currently available only at `1.0 V`; comparing the `0.30..0.80 V` trough requires rerunning or re-exporting Sentaurus with multi-bias plot states and mobility fields.
- Execution of the multi-bias Sentaurus debug-prep step:
  - Local tool availability check: `where.exe sdevice` and `where.exe tdx` do not find Sentaurus tools in the current PATH, so this machine cannot directly rerun/export Sentaurus fields in this session.
  - Added `scripts/prepare_pn2d_sentaurus_multibias_debug.py`.
  - Added regression coverage for generating a Sentaurus multi-bias debug deck.
  - Generated rerun deck:
    - `build/reference_tcad/pn2d_sentaurus2018/sentaurus_debug/pn2d_iv_multibias_debug_sdevice.cmd`
    - `build/reference_tcad/pn2d_sentaurus2018/sentaurus_debug/pn2d_iv_multibias_debug_summary.json`
  - The generated deck keeps the original IV setup and adds:
    - Plot variables `eMobility` and `hMobility`.
    - Quasistationary snapshots: `Plot(FilePrefix="pn2d_iv_multibias" Time=(0;0.25;0.3;0.5;0.8;1) NoOverWrite)`.
  - Once this deck is run in a Sentaurus environment, import each generated `pn2d_iv_multibias*.tdr` snapshot with `build/sentaurus_import --tdr <snapshot.tdr> --export-dir <snapshot_export_dir>`, then rerun the contact-edge transport proxy comparison per bias.
- Inspection of attached `pn2d_iv_multibias_des.tdr`:
  - File read/import succeeded with `build/sentaurus_import.exe`.
  - Exported neutral fields under `build/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias_export/`.
  - HDF5 structure contains only `/collection/geometry_0/state_0`; no multi-state or per-bias collected states were present.
  - ContactExternalVoltage confirms this state is the final `1.0 V` IV point: region 1 `0 V`, region 2 `1 V`.
  - Required final-state fields are present, including `ElectrostaticPotential`, `eDensity`, `hDensity`, `eQuasiFermiPotential`, `hQuasiFermiPotential`, `eCurrentDensity`, `hCurrentDensity`, `TotalCurrentDensity`, `eMobility`, and `hMobility`.
  - Therefore this TDR is sufficient for the `1.0 V` mobility/current-driver branch, but not sufficient for the planned `0.30..0.80 V` trough field-level comparison.
  - Direct `1.0 V` contact-adjacent mobility comparison from this TDR:
    - Cathode electron Sentaurus `0.0727054 m^2/V/s`, Vela `0.0685818 m^2/V/s`, ratio `0.9433`.
    - Anode hole Sentaurus `0.0319098 m^2/V/s`, Vela `0.0321141 m^2/V/s`, ratio `1.0064`.
  - This confirms electron mobility contributes a secondary `~5.7%` current deficit at `1.0 V`; hole mobility is not low. The larger `~13-14%` majority-carrier QF-drop deficit remains the primary observed `1.0 V` transport-driver difference.
- Original Sentaurus IV command update:
  - Modified `reference_tcad/pn2d_sentaurus2018/source/pn2d_iv_sdevice.cmd` so the original rerun can emit independent multi-bias debug TDR snapshots.
  - The command keeps the final `Plot = "pn2d_iv_des.tdr"` and current file output, and adds `eMobility/hMobility` plus `Plot(FilePrefix="pn2d_iv_multibias" Time=(0;0.25;0.3;0.5;0.8;1.0) NoOverWrite)` inside the quasistationary solve.
  - Added regression coverage that locks the original command file to these mobility fields and snapshot times.
- Execution with six numbered Sentaurus multi-bias TDR snapshots:
  - Imported `pn2d_iv_multibias_0000_des.tdr` through `pn2d_iv_multibias_0005_des.tdr` into `build/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias_exports/`.
  - Verified the six states are independent bias points from `ContactExternalVoltage`: `0, 0.25, 0.3, 0.5, 0.8, 1.0 V`.
  - Each imported TDR has the required field set for transport-driver comparison: quasi-Fermi potentials, carrier densities, current densities, and `eMobility/hMobility`.
  - Added reusable multibias contact-edge transport comparison support to `scripts/analyze_pn2d_iv_transport_shape.py`, covered by regression test.
  - Generated reports:
    - `contact_edge_transport_proxy_compare_multibias.csv/json`
    - `contact_edge_transport_proxy_multibias.csv/json`
    - `contact_edge_transport_proxy_multibias_edges.csv`
    - `multibias_transport_driver_key_ratios.csv`
    - `multibias_qf_node_delta.csv`
    - `multibias_weighted_contact_driver.csv`
    - `multibias_weighted_contact_driver_all_carriers.csv`
    - `vela_edge_current_effective_mu_key_bias.csv`
  - Key weighted contact-driver ratios on the Cathode, including both electron and hole branches, track the IV trough:
    - `0.30 V`: weighted driver ratio `0.598`, IV ratio `0.642`.
    - `0.50 V`: weighted driver ratio `0.467`, IV ratio `0.411`.
    - `0.80 V`: weighted driver ratio `0.606`, IV ratio `0.557`.
    - `1.00 V`: weighted driver ratio `0.879`, IV ratio `0.827`.
  - Vela edge-current consistency check rules out a current-extraction or SG coefficient bug for the majority branch:
    - Cathode electron `effective_mu_from_edge_current / mobility` is approximately `1.0000007` at `0.25 V`, `1.0000003` at `0.30 V`, `1.0000006` at `0.50 V`, `0.999986` at `0.80 V`, and `0.996079` at `1.00 V`.
  - Node-level QF split localizes the high-bias difference to the first interior nodes, not the contact boundary values:
    - Cathode electron contact QF delta stays near `0 V`; first interior node delta is `-0.124 mV` at `0.80 V` and `-0.578 mV` at `1.00 V`.
    - Anode hole contact QF delta stays near `0 V`; first interior node delta is `+0.123 mV` at `0.80 V` and `+0.608 mV` at `1.00 V`.
    - Contact carrier densities are aligned; first interior majority density is only about `0.5%` low at `0.80 V` and `1.5-1.6%` low at `1.00 V`.
  - Updated root-cause branch: the dominant remaining difference is a contact-adjacent nonlinear state/QF-gradient deficit in Vela. The current evidence does not support mobility, terminal current extraction, contact QF boundary value, or SG edge-current conversion as the primary cause.
- Execution of the contact-adjacent continuity-balance diagnostic:
  - Added opt-in C++ diagnostic config `sweep.diagnostics.continuity_balance`, with regression coverage in `tests/test_dc_sweep.cpp`.
  - Generated real pn2d IV diagnostic output with VTK disabled:
    - Config: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/simulation_iv_1v_fixed_probe_continuity_fast.json`
    - Raw rows: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed_continuity_fast/iv_1v_fixed_probe_continuity_balance.csv`
    - Key-bias aggregate: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/vela_continuity_balance_key_bias.csv`
  - The sweep converged for 21 points from `0` to `1.0 V`.
  - At key forward biases, contact-edge and neighbor-edge continuity fluxes balance by large cancellation, while recombination is a smaller correction:
    - Cathode electron residual fraction of summed terms: `8.62e-7` at `0.25 V`, `1.58e-8` at `0.30 V`, `1.78e-10` at `0.50 V`, `2.61e-11` at `0.80 V`, and `1.61e-9` at `1.00 V`.
    - Anode hole residual fraction of summed terms: `-4.18e-5` at `0.25 V`, `-1.21e-5` at `0.30 V`, `-6.59e-10` at `0.50 V`, `3.37e-12` at `0.80 V`, and `2.90e-10` at `1.00 V`.
  - The majority-branch contact-adjacent QF drops increase monotonically in Vela (`~1e-9 mV` near `0.25 V`, `~0.188 mV` at `0.80 V`, and `~3.8 mV` at `1.00 V`) and do not show a local numerical residual spike in the IV trough region.
  - Updated branch decision: the local continuity residual itself is not the trough source. The remaining difference is upstream of the residual evaluation: Vela solves a different contact-adjacent quasi-Fermi state than Sentaurus, while satisfying its own discrete continuity equations.
- Execution of solver/path and material-parameter sensitivity probes:
  - Generated solver sensitivity reports under `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/solver_sensitivity/`.
  - Aggregated reports:
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/solver_sensitivity_all_iv_qfdrop.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/multibias_state_field_compare.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/multibias_state_field_signed_compare.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/multibias_srh_recombination_compare.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/multibias_mobility_field_compare.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/sentaurus_inferred_ni_eff_multibias.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/ni_1p45e10_iv_qfdrop_compare.csv`
  - Solver convergence/path probes did not move the IV curve:
    - `reltol=1e-12`, `residual_weights.phin/phip=1e4`, `legacy_node_local`, `warm_start=false`, and `step=0.025 V` all reproduce the baseline ratios to numerical precision.
    - Enabling the old minority-electron contact relaxation makes the IV worse and fails before `1.0 V`; it is not a fix branch.
  - Mobility field parity is not the dominant trough source:
    - Vela Caughey-Thomas electron mobility is about `0.976x` Sentaurus on average; hole mobility is about `1.021x`.
    - The ratios are bias-independent from `0` to `1.0 V`, so mobility cannot explain the bias-local trough shape.
  - SRH recombination is physically mismatched but not a single-parameter IV fix:
    - Recomputed Vela SRH from the solved VTK state is about `35x..59x` larger than Sentaurus over `0.25..1.0 V`.
    - Prior lifetime/recombination scans showed that disabling SRH or sweeping equal lifetimes does not fix the trough without destroying the good low-bias point, so SRH remains a secondary model-parity issue rather than the selected IV-shape fix.
  - Multi-bias state comparison exposes the strongest root-cause candidate:
    - At `0.25/0.30/0.50 V`, global QF differences are small (`~1e-5..1e-4 V`) while carrier densities are systematically low in Vela: mean density ratio is about `0.722..0.725`.
    - Directly inferred Sentaurus effective intrinsic density from `n = ni_eff exp((psi - phin)/Vt)` and `p = ni_eff exp((phip - psi)/Vt)` is stable across all six TDRs: median `ni_eff ~= 1.65562e10 cm^-3`.
    - `scripts/analyze_pn2d_iv_transport_shape.py` now emits `sentaurus_inferred_ni_eff_multibias.csv/json` from exported multibias Sentaurus fields, with regression coverage that reconstructs a known synthetic `ni_eff`.
    - Vela default silicon uses `ni = 1.0e10 cm^-3`; with the current Slotboom factor near `1e17 cm^-3`, Vela effective `ni` is only about `1.13e10 cm^-3`.
    - A one-variable material override `Si ni = 1.45e10 cm^-3` strongly improves the trough/high-bias ratios (`0.50 V: 0.411 -> 0.858`, `0.80 V: 0.557 -> 0.918`, `1.00 V: 0.827 -> 0.966`) but overdrives low bias (`0.25 V: 0.998 -> 1.718`, `0.30 V: 0.642 -> 1.198`).
    - A direct promotion attempt of `Si ni = 1.45e10 cm^-3` into the pn2d_sentaurus2018 IV reference deck regenerated and solved successfully, but failed the existing low-bias comparison gate (`0.2..0.3 V` orders-of-magnitude delta `0.506 > 0.4`). The baseline regenerated deck still passes the same gate (`0.298 < 0.4`). Therefore the material-`ni` parity fix is confirmed as a major root-cause branch but is not sufficient to promote without first fixing the low-bias contact-adjacent QF-drop residual.
  - Updated root-cause branch: the IV mismatch is not caused by current extraction, SG algebra, solver tolerance, sweep path, contact Dirichlet handling, or global mobility. The dominant cause is effective-intrinsic-density/OldSlotboom material parity, with a remaining contact-adjacent QF-drop shape mismatch that was previously masked at `0.25 V` by density/QF-drop error cancellation.
- Execution of Sentaurus-state residual substitution into Vela equations:
  - Added read-only `NewtonSolver::evaluateResidual(state)` and a runner mode `simulation_type = "newton_residual_probe"` that loads external Sentaurus `ElectrostaticPotential/eQuasiFermiPotential/hQuasiFermiPotential` CSV fields and writes node-wise Vela residuals.
  - Real pn2d reports:
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/sentaurus_state_residual_probe/sentaurus_state_residual_probe_status.json`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/sentaurus_state_residual_probe/sentaurus_state_residual_block_compare.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/sentaurus_state_residual_probe/sentaurus_state_residual_top_nodes.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/sentaurus_state_residual_probe/top_continuity_node_inferred_ni_eff.csv`
  - Substituting Sentaurus states into Vela with baseline material gives a large Poisson residual (`~54.3..54.6`) across all biases; using `Si ni = 1.45e10 cm^-3` reduces that Poisson residual to `~1.63`, about `3.0%` of baseline.
  - At the top continuity residual nodes, baseline Vela effective `ni` is consistently about `0.683x` Sentaurus inferred `ni_eff`; the `ni = 1.45e10` material override raises this to about `0.990x`.
  - The phin/phip continuity residuals are much smaller than the Poisson residual at low bias, but scale almost exactly by `1.45x` under the `ni = 1.45e10` override. This explains why material parity fixes the global electrostatic/density state and mid/high-bias trough, while low-bias contact/junction QF-drop current can overdrive: the same tiny QF drops are multiplied by a larger, now Sentaurus-like intrinsic-density factor.
  - The top continuity residuals are not on contact Dirichlet nodes; they concentrate near the junction/contact-adjacent mesh line around `x ~= 1.0 um` and low `y` or midline nodes. This keeps the remaining branch on contact/junction-adjacent QF state shape rather than terminal extraction or boundary value pinning.
- Execution of Sentaurus CSV precision audit:
  - Added a regression test proving that `sentaurus_import` must preserve `~1e-12 V` scalar-field differences in exported CSV files.
  - Updated `SentaurusTdrReader::exportNeutral` field, node, and doping CSV writers to use double round-trip precision instead of the default stream precision.
  - Re-exported `pn2d_iv_multibias_0000_des.tdr` through `pn2d_iv_multibias_0005_des.tdr` and regenerated `analyze_pn2d_iv_transport_shape.py` reports.
  - The high-precision re-export removes the risk that low-bias QF-drop diagnostics are silently rounded away. After re-export, the Cathode/electron contact-edge QF-drop ratios remain: `0.25 V = 1.459`, `0.30 V = 0.867`, `0.50 V = 0.477`, `0.80 V = 0.601`, `1.00 V = 0.868`.
  - Added `sentaurus_contact_edge_inferred_ni_eff_multibias.csv/json`. It shows the Sentaurus contact-edge majority-carrier effective `ni` is flat at `~1.65562e10 cm^-3` across bias and contact/interior locations. The `Si ni = 1.45e10 cm^-3` Vela override corresponds to about `0.9905x` of this local value, so the remaining low-bias overdrive is not explained by contact-edge material `ni` being too high.
  - Updated branch decision: material/OldSlotboom amplitude parity is now well supported, but the remaining unresolved discrepancy is the bias-dependent contact-adjacent QF-drop response. The next targeted diagnostic should compare the local continuity equation coefficients/neighbor-edge cancellation for Vela baseline versus `ni=1.45e10` at `0.25/0.30/0.50 V`, using high-precision Sentaurus QF fields only as a reference state.
- Execution of baseline versus `ni=1.45e10` continuity-balance comparison:
  - Generated `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/solver_sensitivity/ni_1p45e10_continuity/iv_ni_1p45e10_continuity_balance.csv`.
  - Generated aggregate report `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/continuity_balance_baseline_vs_ni145_key_bias.csv`.
  - The `ni=1.45e10` run converged for all 21 sweep points from `0` to `1 V`.
  - At majority contact-adjacent rows, the material override scales QF drop, contact-edge flux, and neighbor-edge flux together; carrier density is unchanged at low/mid bias:
    - `0.25 V`: Cathode/electron QF and both flux terms grow `~1.666x`; Anode/hole grows `~1.580x`; density ratio `~1.0`.
    - `0.30 V`: Cathode/electron grows `~1.807x`; Anode/hole grows `~1.692x`; density ratio `~1.0`.
    - `0.50 V`: Cathode/electron grows `~2.085x`; Anode/hole grows `~2.070x`; density ratio `~1.0`.
    - `0.80 V`: both majority branches grow `~1.65x`; density ratio `~1.005`.
    - `1.00 V`: both majority branches grow `~1.168x`; density ratio `~1.016`.
  - Interpretation: the low/mid-bias overdrive after matching Sentaurus `ni_eff` is not caused by a single continuity term losing balance. Vela's discrete continuity equation remains internally balanced, but its contact-adjacent QF state responds too strongly to the intrinsic-density/material change below about `0.8 V`. The next fix branch should target the equation/model terms that set the local QF state, not contact extraction or SG coefficient algebra.
- Execution of terminal component and spatial QF-response split:
  - Generated:
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/sentaurus_vela_terminal_component_compare.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/sentaurus_vela_terminal_component_compare_loginterp.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/qf_response_baseline_ni145_vs_sentaurus_lowbias.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/qf_response_baseline_ni145_vs_sentaurus_lowbias_summary.csv`
  - Component split shows the baseline low/mid-bias hole current at the current contact is systematically low. With log-interpolated Sentaurus component currents, the Vela/Sentaurus hole ratios are about `0.47` at `0.25/0.30/0.50 V`, increasing to `0.62` at `0.80 V` and `0.858` at `1.00 V`. The `ni=1.45e10` run raises the hole branch to about `0.986x` at `0.25/0.30/0.50 V`.
  - This hole-current behavior matches the contact-edge minority-density ratio: baseline Cathode/hole density is about `0.465x` Sentaurus at `0.25/0.30/0.50 V`, which is `0.683^2` from the baseline effective-`ni` ratio. Therefore the low hole branch is a direct `ni_eff^2` minority-carrier effect.
  - The majority electron branch is different: density is already aligned at the n contact, so the mismatch follows the QF drop. Baseline Cathode/electron contact-edge QF-drop ratios are `1.459` at `0.25 V`, `0.867` at `0.30 V`, and `0.477` at `0.50 V`. With `ni=1.45e10`, those become `2.432`, `1.567`, and `0.995`, respectively.
  - Spatial QF-response report confirms the electron response is not just a terminal extraction artifact:
    - Cathode/electron contact strip (`x` within `0.1 um` of Cathode): `ni=1.45e10` gives Sentaurus ratios `2.432`, `1.567`, `0.995` at `0.25/0.30/0.50 V`.
    - Junction strip (`0.9 <= x <= 1.1 um`): `ni=1.45e10` remains high at `3.609`, `2.578`, `1.504` for the same biases.
  - Updated root-cause split:
    - Hole-current deficit: explained by Vela baseline effective intrinsic density being low; minority density scales approximately with `ni_eff^2`.
    - Electron-current trough/overdrive: not density, mobility, terminal extraction, contact BC value, or SG coefficient; it is the bias-dependent contact-to-junction electron QF-drop shape. Matching `ni_eff` corrects `0.5 V` and high-bias amplitude but over-amplifies the low-bias electron QF drop.

## Files and Responsibilities

- Create: `scripts/analyze_pn2d_iv_transport_shape.py`
  - Reads existing fixed IV VTK files, Sentaurus IV field exports, terminal balance CSV, and contact edge CSV.
  - Emits per-bias transport comparison CSV/JSON and plots.
- Modify: `tests/regression/test_reference_tcad_tools.py`
  - Adds help/smoke tests and pure-function tests for the new diagnostic script.
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_probe_*.vtk`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_probe.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_probe_terminal_balance.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_probe_contact_edges.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/eCurrentDensity_region0.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/hCurrentDensity_region0.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/TotalCurrentDensity_region0.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/ElectricField_region0.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/srhRecombination_region0.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/`
- Candidate modify only after root cause is localized: `src/physics/MobilityModel.cpp`, `src/solver/DriftDiffusionAssembler.cpp`, `src/post/ContactCurrent.cpp`, or `src/simulation/DCSweep.cpp`.

---

### Task 1: Add a Read-Only Transport Shape Diagnostic Script

**Files:**
- Create: `scripts/analyze_pn2d_iv_transport_shape.py`
- Modify: `tests/regression/test_reference_tcad_tools.py`

- [x] **Step 1: Add a help smoke test**

Add this test to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_analyze_pn2d_iv_transport_shape_help(self) -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / "analyze_pn2d_iv_transport_shape.py"),
            "--help",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    self.assertEqual(result.returncode, 0)
    self.assertIn("--reference-root", result.stdout)
    self.assertIn("--biases", result.stdout)
    self.assertIn("--out-dir", result.stdout)
```

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_analyze_pn2d_iv_transport_shape_help
```

Expected: FAIL because the script does not exist yet.

- [x] **Step 2: Create the script CLI and VTK parser**

Create `scripts/analyze_pn2d_iv_transport_shape.py` with:

```python
#!/usr/bin/env python3
"""Compare pn2d IV transport-shape drivers between Sentaurus and Vela."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, default=Path("build/reference_tcad/pn2d_sentaurus2018"))
    parser.add_argument("--biases", default="0.25,0.3,0.5,0.8,1.0")
    parser.add_argument("--out-dir", type=Path, default=Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape"))
    return parser.parse_args()


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    fields: dict[str, list[float]] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {"SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA"}:
                    break
                values.extend(float(v) for v in next_parts)
                i += 1
            fields[name] = values
            continue
        i += 1
    return fields


def abs_stats(values: list[float]) -> dict[str, float]:
    clean = sorted(abs(v) for v in values if math.isfinite(v))
    if not clean:
        return {"points": 0}
    def pct(p: float) -> float:
        idx = (len(clean) - 1) * p / 100.0
        lo = math.floor(idx)
        hi = math.ceil(idx)
        if lo == hi:
            return clean[lo]
        return clean[lo] * (hi - idx) + clean[hi] * (idx - lo)
    return {
        "points": len(clean),
        "mean_abs": sum(clean) / len(clean),
        "median_abs": pct(50),
        "p95_abs": pct(95),
        "max_abs": clean[-1],
    }
```

Run the help test again.

Expected: PASS.

### Task 2: Compare Sentaurus Current-Density Field Magnitudes by Bias

**Files:**
- Modify: `scripts/analyze_pn2d_iv_transport_shape.py`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/current_density_summary.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/current_density_summary.json`

- [x] **Step 1: Add Sentaurus field readers**

Implement a CSV field loader that reads Sentaurus field rows with coordinate columns and magnitude/value columns. If the exact field schema differs, inspect one header and map the numeric data columns explicitly:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import csv
from pathlib import Path
p = Path("build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/eCurrentDensity_region0.csv")
print(next(csv.reader(p.open())))
PY
```

Expected: identify the current-density column names before computing metrics.

- [x] **Step 2: Emit per-bias Sentaurus field magnitude stats**

For each target bias `0.25,0.30,0.50,0.80,1.00 V`, write one row per field:

```text
bias_V,field,points,mean_abs,median_abs,p95_abs,max_abs
```

Fields:

```text
eCurrentDensity
hCurrentDensity
TotalCurrentDensity
ElectricField
srhRecombination
```

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts/analyze_pn2d_iv_transport_shape.py --biases 0.25,0.3,0.5,0.8,1.0
```

Expected: CSV and JSON are generated without changing simulation outputs.

### Task 3: Add Vela-Side Transport Proxies from Existing VTK Fields

**Files:**
- Modify: `scripts/analyze_pn2d_iv_transport_shape.py`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/vela_state_proxy_summary.csv`

- [x] **Step 1: Extract Vela scalar proxy stats**

From each VTK file, compute:

```text
Potential range
ElectronQuasiFermi range
HoleQuasiFermi range
Electrons cm^-3 mean/p95/max
Holes cm^-3 mean/p95/max
NetDoping cm^-3 mean/p95/max
```

Note: Vela VTK carrier densities are in `m^-3`; convert to `cm^-3` by multiplying by `1e-6` before comparing with Sentaurus.

- [x] **Step 2: Compute quasi-Fermi drop proxies**

For each bias, compute:

```text
electron_qf_span_V = max(ElectronQuasiFermi) - min(ElectronQuasiFermi)
hole_qf_span_V = max(HoleQuasiFermi) - min(HoleQuasiFermi)
potential_span_V = max(Potential) - min(Potential)
```

Expected: if current ratio dips while QF span remains smooth, the error likely sits in mobility/flux conversion; if QF span dips similarly, the nonlinear solve state is the source.

### Task 4: Compare Contact Edge Current Concentration Across Bias

**Files:**
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_probe_contact_edges.csv`
- Modify: `scripts/analyze_pn2d_iv_transport_shape.py`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/contact_edge_concentration.csv`

- [x] **Step 1: Aggregate edge concentration**

For each contact and bias, compute:

```text
bias_V
contact
edge_count
sum_current_A_per_um
abs_sum_current_A_per_um
max_abs_edge_current_A_per_um
top3_abs_fraction
edge_sign_changes
```

Expected: if `top3_abs_fraction` spikes in the `0.3..0.8 V` trough, debug contact-edge geometry and boundary integration.

- [x] **Step 2: Compare against terminal balance**

Join with `iv_1v_fixed_probe_terminal_balance.csv` and verify:

```text
abs(edge_sum - terminal_total) < 1e-18 A/um
```

Expected: keep current extraction ruled out unless this gate fails.

### Task 5: Add C++ Transport Diagnostics Only if Python Proxies Are Insufficient

**Files:**
- Modify if needed: `src/simulation/DCSweep.cpp`
- Modify if needed: `src/solver/DriftDiffusionAssembler.cpp`
- Test if modified: `tests/test_dc_sweep.cpp`

- [x] **Step 1: Write a failing test for opt-in transport columns**

If Python diagnostics cannot isolate the layer, add a test that enables:

```json
"sweep": {
  "diagnostics": {
    "transport": { "enabled": true }
  }
}
```

Assert output CSV contains:

```text
mean_electron_mobility_m2_V_s
mean_hole_mobility_m2_V_s
min_electron_mobility_m2_V_s
min_hole_mobility_m2_V_s
max_electric_field_V_per_cm
mean_electron_qf_gradient_V_per_cm
mean_hole_qf_gradient_V_per_cm
```

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
ctest --test-dir build --output-on-failure -R DCSweep
```

Expected before implementation: FAIL due to missing columns.

- [x] **Step 2: Implement opt-in diagnostics only**

Keep the default CSV schema unchanged unless `sweep.diagnostics.transport.enabled=true`.

Expected: existing reference outputs remain stable by default.

### Task 6: Decision Gates for the Next Fix

**Files:**
- Candidate modify after evidence:
  - `src/physics/MobilityModel.cpp`
  - `src/solver/DriftDiffusionAssembler.cpp`
  - `src/post/ContactCurrent.cpp`
  - `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`

- [x] **Step 1: Select exactly one root-cause branch**

Use these gates:

```text
If current-density field ratio dips in Sentaurus but Vela terminal current dips more, debug ContactCurrent/contact-edge integration.
If QF span and carrier density match but Vela current remains low, debug SG flux assembly and mobility/Einstein relation.
If QF span is already low in Vela in the same bias range, debug nonlinear solve state and boundary conditions.
If SRH recombination changes sharply in Sentaurus around the trough and Vela lacks the same trend, debug recombination/effective-ni implementation.
If all proxies are smooth while only terminal current dips, debug unit convention or width/current-density conversion one layer deeper.
```

- Execution of `ni=1.45e10` SRH lifetime sweep after the contact-edge QF overdrive was isolated:
  - Ran four single-variable Vela sweeps with the already-aligned material `ni` and unchanged mobility/mesh/contact settings:
    - `taun=taup=1e-6 s`
    - `taun=taup=2e-6 s`
    - `taun=taup=4e-6 s`
    - `taun=taup=6e-6 s`
  - Generated:
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/ni145_srh_lifetime_sweep_compare.csv`
    - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/ni145_srh_lifetime_sweep_summary.csv`
  - The comparison uses the same linear-current interpolation convention as the existing IV gate and `ni145_no_srh_iv_qfdrop_compare.csv`.
  - Key Vela/Sentaurus total-current ratios:
    - Baseline context: `0.25V 0.998`, `0.30V 0.642`, `0.50V 0.411`, `0.80V 0.557`, `1.00V 0.827`.
    - `ni=1.45e10, tau=1e-7`: `0.25V 1.718`, `0.30V 1.198`, `0.50V 0.858`, `0.80V 0.918`, `1.00V 0.966`.
    - `ni=1.45e10, tau=1e-6`: `0.25V 0.922`, `0.30V 0.877`, `0.50V 0.848`, `0.80V 0.918`, `1.00V 0.967`.
    - `ni=1.45e10, tau=2e-6`: `0.25V 0.877`, `0.30V 0.859`, `0.50V 0.847`, `0.80V 0.918`, `1.00V 0.967`.
    - `ni=1.45e10, tau=4e-6`: `0.25V 0.855`, `0.30V 0.850`, `0.50V 0.847`, `0.80V 0.918`, `1.00V 0.967`.
    - `ni=1.45e10, tau=6e-6`: `0.25V 0.848`, `0.30V 0.847`, `0.50V 0.847`, `0.80V 0.918`, `1.00V 0.967`.
    - `ni=1.45e10, recombination=none`: `0.25V 0.833`, `0.30V 0.841`, `0.50V 0.846`, `0.80V 0.918`, `1.00V 0.967`.
  - Cathode electron branch is the sensitive branch:
    - At `0.25V`, `ni=1.45e10, tau=1e-7` has electron current `1.650e-14 A/um` and mean Cathode electron QF drop `9.278e-10 mV`.
    - At `0.25V`, `tau=1e-6` lowers these to `7.456e-15 A/um` and `4.192e-10 mV`.
    - At `0.25V`, `recombination=none` lowers these further to `6.451e-15 A/um` and `3.626e-10 mV`.
  - Branch decision update: the earlier material-`ni` mismatch is still the dominant amplitude/root electrostatics mismatch, but once material `ni_eff` is aligned, default Vela SRH (`taun=taup=1e-7 s`) overdrives the low-bias contact-adjacent electron QF response. Increasing equal SRH lifetimes to about `1e-6 s` removes most of that low-bias overdrive while preserving the mid/high-bias improvement. Larger lifetimes approach the no-SRH limit and under-shoot the low-bias current.
  - Follow-up SRH field-magnitude diagnostic:
    - Enabled `solver.diagnostics=true` for three `ni=1.45e10` sweeps and generated `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/ni145_srh_recombination_diagnostic_compare.csv`.
    - `taun=taup=1e-7 s` gives Vela/Sentaurus mean-abs SRH ratios of about `54x..64x` across `0.25..1.00 V`, matching the low-bias electron-QF overdrive symptom.
    - `taun=taup=1e-6 s` gives better IV ratios (`0.922`, `0.877`, `0.848`, `0.918`, `0.967`) but still leaves mean-abs SRH about `5.4x..6.4x` Sentaurus.
    - `taun=taup=6e-6 s` matches SRH field magnitude much better: mean-abs ratio `0.905x..1.06x`, max-abs ratio `1.07x..1.08x`, while IV ratios are `0.848`, `0.847`, `0.847`, `0.918`, `0.967`.
  - Refined root-cause statement: the IV discrepancy is not a single-parameter error. Baseline material/effective-intrinsic-density is too low, which suppresses minority density and the mid/high forward current. After aligning `ni_eff`, the default Vela SRH lifetime is too short by roughly `~60x` relative to the Sentaurus SRH field magnitude, which overdrives the low-bias electron quasi-Fermi drop. A lifetime near `1e-6 s` is an IV-fitting compromise, while a lifetime near `6e-6 s` is closer to SRH field parity.
  - Execution of direct `models.par` parameter-alignment experiment:
    - User supplied `reference_tcad/pn2d_sentaurus2018/source/models.par`, generated by `sdevice -P:Silicon`.
    - Extracted Sentaurus Silicon defaults relevant to this IV deck:
      - `Bandgap`: `Eg0=1.16964 eV`, `dEg0(OldSlotboom)=-1.595e-2 eV`, `alpha=4.73e-4 eV/K`, `beta=636 K`.
      - `OldSlotboom`: `Ebgn=9e-3 eV`, `Nref=1e17 cm^-3`, `C=0.5`.
      - `Scharfetter` SRH: electron `taumax=tau0=1e-5 s`; hole `taumax=tau0=3e-6 s`; `Etrap=0`.
    - Generated a diagnostic Vela run under `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/solver_sensitivity/models_par_alignment/` with:
      - material `Si ni=1.4638914958767616e10 cm^-3`;
      - existing Vela Slotboom parameters, which already match the Sentaurus `OldSlotboom` block;
      - `taun=1e-5 s`, `taup=3e-6 s`;
      - `solver.diagnostics=true` and contact-edge diagnostics enabled.
    - Generated comparison reports:
      - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/models_par_alignment_iv_srh_compare.csv`
      - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/models_par_alignment_iv_srh_summary.csv`
    - Key total-current ratios for `models_par_alignment`: `0.25V 0.864`, `0.30V 0.863`, `0.50V 0.863`, `0.80V 0.929`, `1.00V 0.971`.
    - SRH field magnitude is essentially matched: mean-abs SRH ratio `1.0002`, `1.0002`, `0.9998`, `1.0000`, `1.0010` at `0.25/0.30/0.50/0.80/1.00 V`; max-abs ratio is also about `1.000x`.
    - Branch decision update: `models.par` explains both the Sentaurus effective intrinsic-density level and the SRH field-magnitude mismatch. With those material/SRH parameters aligned, the remaining IV error is a smoother `~0.86x` low/mid-bias current deficit rather than the original trough/overdrive shape. Therefore the next unresolved branch is no longer material `ni` or SRH lifetime; it is most likely the residual transport discretization/mobility-detail layer, especially full Masetti `DopingDependence` versus Vela's simplified Caughey-Thomas-like mobility or the SG/Einstein coefficient implementation.
  - Execution of the mobility/comparison-metric branch:
    - Generated `models_par_alignment_no_field_mobility` by changing only `solver.mobility.model` from `caughey_thomas_field` to `caughey_thomas`. This changed IV by only `~0.2%..0.5%`, so high-field mobility limiting is not the remaining low/mid-bias root cause.
    - Generated `models_par_alignment_constant_sentaurus_mobility` with `mobility=model constant` and material mobility set to the Sentaurus contact-adjacent values: `mun=727.054 cm^2/V/s`, `mup=319.098 cm^2/V/s`.
    - Generated reports:
      - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/models_par_alignment_mobility_branch_compare.csv`
      - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/models_par_alignment_mobility_branch_summary.csv`
      - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/models_par_alignment_mobility_loginterp_component_compare.csv`
      - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/models_par_alignment_mobility_loginterp_component_summary.csv`
    - Against the old linear-current interpolation reference, constant Sentaurus mobility gives ratios `0.25V 0.894`, `0.30V 0.894`, `0.50V 0.893`, `0.80V 0.958`, `1.00V 1.001`.
    - Against log-current interpolation/component-current references, the same run gives:
      - Total component ratio: `0.25V 0.999`, `0.30V 1.000`, `0.50V 1.000`, `0.80V 1.023`, `1.00V 1.001`.
      - Electron component ratio: `0.999`, `1.000`, `1.000`, `1.023`, `1.001`.
      - Hole component ratio: `1.000`, `1.000`, `1.000`, `1.023`, `1.000`.
    - The old linear interpolation reference is higher than the log/current-component reference by `1.117x`, `1.119x`, `1.120x`, `1.068x`, and `1.000x` at `0.25/0.30/0.50/0.80/1.00 V`.
  - Final branch decision update: after `models.par` material/SRH parameters and Sentaurus-equivalent mobility are aligned, Vela matches Sentaurus component currents to about `0.1%` through `0.5 V` and within about `2.3%` at `0.8 V`. The apparent `~0.89x` residual under the older table is primarily an artifact of linear current interpolation on an exponential IV curve, combined with the simplified Vela mobility being about `4.8%` low for the contact-adjacent electron branch. The remaining actionable implementation work is to (1) promote Sentaurus-derived material/SRH defaults/config generation, (2) add a Masetti/Sentaurus `DopingDependence` mobility model or equivalent parameter import, and (3) update IV comparison tools to use log-current interpolation for current magnitude checks.
  - Execution of the first implementation hardening step:
    - Added explicit `log_current` interpolation support to `scripts/compare_reference_curves.py`. The historical default remains `linear`; `log_current` only applies when requested and falls back to linear interpolation across zero or sign-changing adjacent samples.
    - Threaded `comparison.interpolation` through `scripts/sentaurus_import.py` so reference configs can opt into log-current comparison gates.
    - Added `reference_tcad/pn2d_sentaurus2018/source/pn2d_sentaurus2018_iv_materials.json` with the `models.par`-derived Silicon intrinsic density `ni=1.4638914958767616e10 cm^-3`.
    - Updated the pn2d Sentaurus2018 IV simulation config to use that material file, `taun=1e-5 s`, `taup=3e-6 s`, and `comparison.interpolation="log_current"`.
    - Regenerated the pn2d Sentaurus2018 reference workspace with the updated config. The IV comparison gate now reports `status=pass`, `interpolation=log_current`, `orders_of_magnitude=0.0147648`, and `max_relative_error=0.0334257` over the configured `0.2..0.3 V` gate window.
  - Execution of the Masetti implementation branch:
    - Added a C++ Masetti mobility model matching Sentaurus Silicon
      `DopingDependence` Formula 1 defaults from `models.par`.
    - Added regression coverage for the exact `1e17 cm^-3` pn2d doping point:
      electron mobility `0.07270544030120773 m^2/(V s)` and hole mobility
      `0.03190980929489245 m^2/(V s)`.
    - Updated Sentaurus physics import so `DopingDependence` maps to
      `mobility.model = masetti`; `DopingDependence + HighFieldSaturation`
      maps to `masetti_field`.
    - Updated the pn2d Sentaurus2018 IV deck from `caughey_thomas_field` to
      `masetti`, because the cmd physics block contains `DopingDependence` but
      no `HighFieldSaturation`.
    - Regenerated the pn2d Sentaurus2018 reference workspace. The IV comparison
      gate now reports `status=pass`, `interpolation=log_current`,
      `orders_of_magnitude=0.000203179`, and `max_relative_error=0.000467947`
      over the configured `0.2..0.3 V` gate window.
    - Branch decision update: after `models.par` material/SRH alignment,
      log-current comparison, and full Masetti `DopingDependence`, the low-bias
      IV discrepancy is reduced to below `0.05%` in the configured gate window.
      This confirms the prior `~3.3%` residual was primarily the simplified
      Caughey-Thomas mobility shape and the accidental high-field model
      selection, not recombination or material `ni`.

- [ ] **Step 2: Add one failing regression before changing physics**

Deferred to the next implementation plan. This execution selected the branch but did not change physics behavior.

For the selected branch, add a minimal test that reproduces the measured discrepancy. Do not change multiple physics axes in one commit.

### Task 7: Verification

**Files:**
- Test: `tests/regression/test_reference_tcad_tools.py`
- Test if C++ modified: `tests/test_dc_sweep.cpp`

- [x] **Step 1: Run Python regression tests**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools
```

Expected: all tests pass.

- [x] **Step 2: Run focused C++ tests if transport diagnostics were added**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "DCSweep|ContactCurrent|SG|DriftDiffusion"
```

Expected: focused tests pass.

- [x] **Step 3: Regenerate IV comparison plot**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts/analyze_pn2d_iv_transport_shape.py --biases 0.25,0.3,0.5,0.8,1.0
```

Expected: reports under `build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape/` identify one next branch.
