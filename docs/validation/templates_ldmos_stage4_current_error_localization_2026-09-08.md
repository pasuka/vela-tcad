# Templates/LDMOS Stage-4 D5 current-error localization

Date: 2026-09-08 (Asia/Shanghai). Worktree `templates-ldmos-phase-a`, source
base `959c1c9` plus the pre-existing optional step-growth changes.

## Result

The edge-projected quasi-Fermi driving field in high-field mobility is a
major contributor to the Vg=8 current excess. At Vd=4 V, switching only the
existing gradient-discretization option to `transport_cell_vector` and
obtaining a new self-consistent solution reduces the current error from
**12.40768% to 1.52738%**, under all original residual, local-row and KCL gates.
This is local causal evidence, not qualification of a replacement full curve.

At Vd=40 V the fixed-state evidence points to the same transport term, but
an abrupt switch from the old accepted 40 V state fails two bounded
160-update attempts. No candidate 40 V current is accepted or scored.
The existing full D5 verdict remains **fail**; IALMob remains off.

No production C++ source, material parameter, frozen reference, acceptance
threshold or previously accepted curve was changed in this investigation.
Generated probes, local candidate configurations and plots are isolated in
`reference_staging/templates_ldmos_current_error_20260908/`.

## Configuration and controls

- UCRT64 Release, Eigen SparseLU/COLAMD with L2 row/column equilibration.
  Runner SHA256 `9638411ec0d38434f05afc765d573060adde0d406c6818f3428c76777a066391`.
- Exact mesh: 10,241 nodes, 19,782 triangles, 30,022 edges; 5,723 silicon
  nodes. Carrier transport uses the same external AverageBox couples,
  barycentric volumes and material-local Poisson charge. Current units A/um,
  potentials V, stored carrier densities m^-3, mesh coordinates um.
- D5 at 300 K: Fermi/OldSlotboom, SRH/Auger, constant-field HFS and live
  mobility Jacobian derivatives; predictor, IALMob and avalanche off.
- Vela states are from today's independent full Vg=8 sweep at 4 and 40 V.
  Sentaurus states are the saved independent 4 V endpoint and adaptive-on
  40 V endpoint. The HDF5 importer exported the latter successfully; all
  10,241 coordinates match exactly. No new remote simulation was needed.
- The 40 V Vela state is stored in the validated +28 V inverse-translation
  frame (raw source -28, gate -20, drain 12 V). Fixed-state field comparison
  adds 28 V to physical potentials; inactive carrier placeholders stay zero.
  The local 40 V Newton control retains the original raw frame.
- Hash checks cover the runner, mesh, material/doping inputs, external
  couples, states and exported fields. Both endpoint input audits pass.

The production `newton_carrier_term_probe` uses
`carrier_term_probe.solved_equation_terms=true`; SG edge and Poisson term
probes use the same input configuration. A probe's successful return means
diagnostics completed, not that its foreign state is a converged solution.
The 5,674 free silicon nodes exclude only configured electrical electrodes.
The mesh's `th_lat` thermal boundary is not an electrical terminal and is
not removed from that free-node set or included in electrical KCL.

## 1. Terminal evaluation versus internal transport

| Vg=8 endpoint | Vd=4 V | Vd=40 V |
| --- | ---: | ---: |
| Sentaurus reference current, mA/um | 0.2271551512 | 0.4190277906 |
| Original accepted Vela current, mA/um | 0.2553398409 | 0.4869356376 |
| Original current error | +12.40768% | +16.20605% |
| Vela SG drain cut evaluated on the Sentaurus state, mA/um | 0.2272046352 | 0.4191313676 |
| Same-state drain-cut error | +0.02178% | +0.02472% |

The drain reconstruction itself reproduces the reference much more closely
when supplied the reference internal state. This narrows the dominant
discrepancy to the different internal solution, rather than a uniform current
unit/width factor or drain summation error. The foreign-state cuts do not
satisfy strict global continuity and are not substitutes for accepted solves.

Internal graph cuts use all transport edges crossing a specified lateral
coordinate y, with electron and hole conventional-current signs. On the
original accepted Vela state, their absolute current agrees with its drain
current (4 V within about 1.1e-13 relative; translated 40 V within about 6e-10).
On the Sentaurus state the edge-HFS formula instead produces a localized
excess near y=4 um:

| Internal y=4 um cut error relative to Sentaurus terminal current | Vd=4 V | Vd=40 V |
| --- | ---: | ---: |
| Edge-projected HFS | +18.1023% | +26.0545% |
| Cell-vector HFS, same state | -3.2101% | -5.1635% |

The vector candidate reduces the major excess but creates a smaller deficit
there. It does not exactly reproduce the proprietary operator.

## 2. Location and term responsible

Changing only `solver.mobility.high_field_gradient_discretization` from
`edge_projection` to `transport_cell_vector` on the frozen Sentaurus states
gives the following production electron residuals on free silicon nodes.
Norms are Vela's internal scaled diagnostic norms, not SDevice RHS values.

| Diagnostic | Vd=4 V | Vd=40 V |
| --- | ---: | ---: |
| Edge HFS electron L2 | 5.195937 | 20.184197 |
| Vector HFS electron L2 | 1.268857 | 5.461116 |
| Vector / edge L2 | 0.244202 | 0.270564 |
| Squared baseline residual fraction in x<-9.5 um, 3.4<=y<4.5 um | 94.67% | 88.27% |

The dominant region is the curved Si/oxide transition under the gate. The
largest baseline electron row at both voltages is node 3615,
(x,y)=(-9.965914,3.641602) um. Candidate hotspots move to node 3432 at 4 V
and node 3726 at 40 V; residual redistribution is retained in the evidence.
The electron recombination-term L2 is only about 4.1e-13 and 5.3e-13 of
the respective flux-term L2 on these frozen states. This screening evidence
does not support SRH/Auger as the dominant source of the observed hotspot.

Current source interpretation:

- [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp) selects
  an edge QF difference divided by edge length, or the recovered cell-vector
  field, before mobility evaluation.
- [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) applies
  the existing contact electric-field fallback. It takes precedence over
  either bulk QF option, explaining why the fixed-state drain cut is
  identical for those two variants.
- [MobilityModel.cpp](../../src/physics/MobilityModel.cpp) applies
  `mu0/[1+(mu0*F/vsat)^beta]^(1/beta)`. Along-edge projection can omit a
  transverse QF gradient on a non-aligned 2-D mesh and therefore weaken
  velocity saturation. The existing cell-vector option reconstructs the
  magnitude using semiconductor P1 cell gradients and adjacent-cell weights.
- The retrieved T-2022.03-SP2 Silicon parameter file has alpha=0,
  mu0=1417/470.5 cm^2/(V s), beta=1.109/1.213 and
  vsat=1.07e7/8.37e6 cm/s at 300 K, consistent with the configured limiter.
  Geometry-based driving-field recovery remains a distinct choice.

Turning contact fallback off at 4 V worsens free-silicon electron L2 to
16.66363 and gives a 2.5632 mA/um foreign-state drain cut. Removing HFS
entirely worsens L2 to 27.04668. Neither is a viable correction.
The prior low-gate fixed-state HFS experiment also found incomplete
improvement; today's result does not retroactively clear its failed gates.

## 3. Self-consistent local test

At Vg=8, Vd=4 V, the vector candidate starts from the accepted baseline
checkpoint. The first bounded attempt applies 124 updates and is rejected
at line-search stagnation (electron norm approximately 2.09e-11, ceiling
1e-11). A strict same-bias reclose applies another 14 updates and passes:

| Accepted candidate check | Result | Original limit |
| --- | ---: | ---: |
| Psi block | 9.84799e-10 | 5e-8 |
| Electron block | 1.87403e-12 | 1e-11 |
| Hole block | 3.40346e-25 | 3e-10 |
| Local carrier-row maximum | 9.78709e-9; zero violations | 1e-8 |
| KCL/max-terminal current | 2.22022e-14 | 1e-8 |

All 10,241 state records are unique and finite, carrier densities are
nonnegative, and contact voltage differences remain Vg=8/Vd=4 V.
No density recovery or predictor is used. The combined local attempt wall
time is 83.530 s; this is not a voltage-sweep performance benchmark.

The accepted current is **0.2306246631 mA/um**, error **+1.52738%**.
The original excess falls by about 87.7%. This local intervention provides
stronger evidence than a fixed-state residual reduction alone.

The 40 V direct-switch experiment and its strict reclose both exhaust 160
updates. They are rejected; the reclose ends near psi=1.42e-4,
electron=2.52e-3 and hole=1.37e-8, with local-row violations still present.
The reported zero current in the failed curve row is a failure placeholder,
not a physical result. No high-voltage candidate score is inferred.

## 4. Field mapping and remaining uncertainty

Across all silicon nodes, the original accepted-state potential difference
has P95/max 0.2182/0.2516 V at Vd=4 V and 3.1187/4.1967 V at Vd=40 V.
Electrostatic and electron-QF changes largely track each other; the maximum
change of psi-phin is 0.09594 V and 0.28082 V, respectively. These are
spatial solution differences, not an arbitrary fitted gauge shift.
Contact QFs agree to roundoff and contact psi discrepancies are sub-mV.

Reconstructing density from the Sentaurus potentials with Vela statistics
has maximum multiplicative discrepancy about 0.512% for electrons and
0.605% for holes. The associated QF mapping discrepancies are sub-mV.
This is smaller than the current error, but it cannot be dismissed globally:
small density errors can dominate Poisson cancellation in heavily doped
regions. At 4 V, a diagnostic remapping to density-consistent QFs reduces
free-silicon Poisson L2 from 1243.0 to 85.94 but generates large carrier
residuals near high-conductivity contacts. That artificial state is not an
accepted solution and does not establish a correction to the model.

An additional import artifact was isolated: 4,065 non-transport nodes in
each archived Sentaurus state carry nonzero carrier QF placeholders. Vela's
inactive carrier rows require zero, so these values dominate the full raw
carrier norm (about 9865 at 4 V) despite carrying no physical transport.
Setting only those inactive phin/phip entries to zero leaves all silicon
fields unchanged and removes the artificial gauge rows. The 4 V total
electron norm then equals the physical 5.19594 baseline, or 1.26886 for
vector HFS. Poisson remains unsatisfied. This corrects the interpretation
of the earlier failed imported-seed control; it does not prove that its
entire failure was due to the placeholders. No importer production code was
changed, and no normalized seed was promoted to an accepted state.

![Field differences and internal cuts](../../reference_staging/templates_ldmos_current_error_20260908/current_error_localization.png)

The spatial plots show silicon nodes only, with lateral coordinate y and
depth coordinate x+10, cropped to depth<1.5 um. Each color range is symmetric
and clipped at the 98th percentile of the plotted absolute difference;
reported maxima above are computed from all silicon nodes. Density log
ratios use the documented 1 m^-3 floor. The internal-cut panels are frozen
operator evaluations, not additional simulated Id-Vd curves.

## Follow-up boundary

Next, continue the qualified vector-HFS low-bias state gradually to the
remaining exact Vg=8 reference points, then verify a fresh zero-to-40 V path
and Vg=4 non-regression before considering a production default change.
Investigate the remaining interface stencil and statistics/Poisson mismatch
separately. Normalize inactive carrier placeholders in future diagnostic
seed preparation with explicit invariance checks. Full two-gate D5 and the
40 V ratio must pass their unchanged gates before IALMob work.

## Evidence and verification

- [Machine-readable findings and candidate gate audit](../../reference_staging/templates_ldmos_current_error_20260908/findings.json)
- [4 V full fixed-state evidence](../../reference_staging/templates_ldmos_current_error_20260908/summary.json)
- [40 V full fixed-state evidence](../../reference_staging/templates_ldmos_current_error_20260908/vd40/summary.json)
- [4 V accepted candidate](../../reference_staging/templates_ldmos_current_error_20260908/vector_reclose_vd4_retry/curve.csv)
- [40 V rejected reclose](../../reference_staging/templates_ldmos_current_error_20260908/vector_reclose_vd40_retry/result.json)
- [4 V inactive-placeholder control](../../reference_staging/templates_ldmos_current_error_20260908/inactive_placeholder_control/summary.json)
- [40 V inactive-placeholder control](../../reference_staging/templates_ldmos_current_error_20260908/vd40/inactive_placeholder_control/summary.json)
- [Probe driver](../../reference_staging/templates_ldmos_current_error_20260908/analyze.py),
  [local-solve driver](../../reference_staging/templates_ldmos_current_error_20260908/reclose_vector.py),
  [placeholder control](../../reference_staging/templates_ldmos_current_error_20260908/normalize_inactive.py),
  [audit and plot generator](../../reference_staging/templates_ldmos_current_error_20260908/summarize.py).

Thirty-three fixed-state probe executions and four bounded coupled attempts
were completed. The latter contain one accepted final candidate (4 V
strict reclose) and three rejected attempts. Audit scripts check source/input
hashes, candidate original gates, state integrity, bias and model settings.
All output files and failures are retained. Scientific figures were inspected.
No solver rebuild or CTest rerun was needed for this analysis-only change.

Reproduce in a fresh evidence directory or preserve the existing probe cache:

```powershell
D:\msys64\ucrt64\bin\python.exe reference_staging/templates_ldmos_current_error_20260908/analyze.py
D:\msys64\ucrt64\bin\python.exe reference_staging/templates_ldmos_current_error_20260908/analyze.py --density-mapping-control
D:\msys64\ucrt64\bin\python.exe reference_staging/templates_ldmos_current_error_20260908/analyze.py --bias 40
D:\msys64\ucrt64\bin\python.exe reference_staging/templates_ldmos_current_error_20260908/normalize_inactive.py
D:\msys64\ucrt64\bin\python.exe reference_staging/templates_ldmos_current_error_20260908/summarize.py
```

The 40 V command requires the exported endpoint TDR under `vd40/export`.
The local-solve driver takes `--bias 4|40` and optional `--strict-reclose`;
it refuses to overwrite completed attempt results.
