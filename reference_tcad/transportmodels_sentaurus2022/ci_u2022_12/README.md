# TransportModels DD and DG regressions

## Id–Vd

`reference_u2022_transportmodels_idvd` runs `pp12_des.cmd` through the built
`sd-sim` with the Vela `sentaurus_u2022_12` constants profile and single-threaded
oneMKL PARDISO. Set `VELA_MKL_RT_LIBRARY` to the absolute path of the MKL runtime.

```sh
cmake -S . -B build-u2022 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DVELA_PHYSICAL_CONSTANTS_PROFILE=sentaurus_u2022_12 \
  -DVELA_ENABLE_MKL_PARDISO=ON -DVELA_ENABLE_TRANSPORTMODELS_IDVD_CI=ON
cmake --build build-u2022 --target sd-sim vela_example_runner
ctest --test-dir build-u2022 -R '^reference_u2022_transportmodels_idvd$' --no-tests=error --output-on-failure
```

The fixture stages GRD/DAT from the archived 3,315 vertices and 6,456 triangles
in `../vela/mesh.json` and the separately exported donor and acceptor
concentrations in `../vela/doping.csv`. The [artifact manifest](../artifact_manifest.json)
identifies the original reflected `n1_msh.tdr`. The test uses the tracked neutral
mesh and does not require a TDR reader or a Sentaurus installation.

The [archived CMD](../source/pp12_des.cmd) supplies the model options and
`Solve` sequence. The test deck uses the GRD grid path, explicit isothermal
300 K, and final potential output. It omits unused thermal boundaries and a
file deletion command. The [archived silicon parameter file](../source/pp12_des.par)
is staged with the PolySi parameters extracted from `sdevice -P:PolySilicon`
in Sentaurus T-2022.03-SP2. The [fixture manifest](manifest.json) pins both the
original export and the staged parameter excerpt. At 300 K the imported
coefficients give Nc = 2.858252e19 cm⁻³, Nv = 3.104629e19 cm⁻³ and an intrinsic
density of 1.463881e10 cm⁻³ before local doping-dependent narrowing. The
sidecar requests unit-scaled electrostatic potential and carrier densities,
with Newton backtracking.

At Vg = 1 V, the deck initializes with up to 100 Poisson iterations followed
by up to 20 coupled iterations. The drain sweep from 0 to 2 V uses
normalized time steps of 0.001, 0.00001 and 1 for initial, minimum and maximum,
respectively. Increment and Decrement default to 2. Extrapolate is enabled.
Vela retains the preceding concentration at nodes below 1e-50 cm⁻³ to
prevent extrapolation near the numerical positivity floor. A local Newton
trial that would turn a concentration nonpositive retains half of its
preceding positive value.
The next step restores the nominal proposal after a successful output-only
clip, with 21 requested output biases from 0 to 2 V in 0.1 V increments.

`sentaurus_reference.csv` is the DD Id–Vd subset of the native currents
recorded in the [historical comparison](../../../docs/validation/transportmodels_dd_dg_srh_corrected_cold_regression_2026-08-23.json).
Its source is Sentaurus T-2022.03-SP2 `IdVds_n12_des.plt`. The source
artifact SHA-256 is recorded in the archive manifest. Bias is in volts and
current is in A/µm.

The test compares the requested biases directly and applies the existing
TransportModels contract of 5% maximum relative drain-current error over
0.1–2 V and 3% at 2 V. It records the near-zero 0 V current and four-terminal
current sum for inspection. The results appear in `result.json` and
`comparison.csv` under the run's CTest report directory.

This is an opt-in qualification test. With the previous PolySilicon parameter
fit, C560 Linux and T470p Windows Release runs passed the 0.1, 0.2 and 0.3 V
current comparisons with relative errors of 1.10%, 1.74% and 2.35%. Those runs
reached the minimum-step limit near 0.30–0.32 V during coupled backtracking.
In a Windows Release comparison on the same mesh, the native coefficients
reduced those three errors to 0.179%, 0.093% and 0.301%, respectively. The
coupled sweep still reduced its step below the minimum near 0.30 V, leaving
the 2 V endpoint unqualified.

With the pinned original mesh and current sidecar, the Release qualification
completed all 21 points. The sidecar selects cell-material permittivity,
Delaunay transfer, signed Poisson charge support, and element-box transport.
It freezes high-field mobility feedback in the Jacobian. The vertex-sampled
Masetti path and effective-potential Fermi–Dirac SG gave 0.034% maximum drain
current error, including the 2 V endpoint, with a maximum terminal current
sum magnitude of 9.10e-13 A/µm.

## Id–Vg

Enable `VELA_ENABLE_TRANSPORTMODELS_IDVG_CI=ON` with the same Release,
physical-constants and oneMKL settings as Id–Vd. The CTest name is
`reference_u2022_transportmodels_idvg`. It executes the staged
`pp6_des.cmd` through `sd-sim` on the same 3,315-node grid and appends the
same PolySi excerpt to the original `pp6_des.par`.

The standalone CMD preserves the original isothermal physics and Solve
sequence. It replaces the TDR grid path with GRD, explicitly sets 300 K,
omits the unused thermal boundaries and output deletion command, and limits
Plot fields to potential. The run starts at `Vg=-1 V, Vd=0 V`, raises the
drain to 1.1 V, then sweeps the gate to 2.2 V with the CMD step bounds and
21 exact CurrentPlot biases. Each sweep consumes the preceding stage state.
The complete sweep uses potential and carrier concentrations. A short
independent sweep from the same equilibrium and drain ramp covers the
first three gate biases through -0.68 V. It uses local contact-basin
quasi-Fermi references to resolve the leakage currents and preserves the
original gate sweep's absolute initial and minimum steps and requested
output locations. The maximum step is capped at the short sweep's span.
The deep-off gate uses the directly measured four terminal currents from
the short sweep. The report also retains the full sweep's drain current at
these biases.

The 21-point oracle is the native DD curve recorded in the
[historical comparison](../../../docs/validation/transportmodels_dd_dg_srh_corrected_cold_regression_2026-08-23.json).
The manifest pins the original Sentaurus T-2022.03-SP2
`IdVgs_n6_des.plt` SHA-256, source CMD and PAR, and the extracted CSV.

The [existing contract](../vela/contracts/transportmodels_dd_dg_sentaurus2022_v1.json)
sets a maximum error of 0.15 dex from deep off through transition and 10%
relative error in the on region. The three deep-off points also require
`Id / abs(sum(Iterminal)) >= 10`. The test reads the four terminal currents
and requested biases directly from the published PLT and writes the
per-bias comparison and failed bias list under `reference-u2022-reports`.

## DG Id–Vd

Enable `VELA_ENABLE_TRANSPORTMODELS_DG_IDVD_CI=ON` to register
`reference_u2022_transportmodels_dg_idvd`. Use the same Release,
physical-constants and oneMKL settings as DD Id–Vd, with HDF5 state
support and Python h5py. The test stages the original `pp13_des.par`
and the shared PolySi excerpt on the same frozen grid.

The standalone CMD keeps `eQuantumPotential` in Physics and in the Solve
statements. The sidecar selects the `sentaurus_box` discretization from the
DD/DG contract. Vela initializes at
`Vg=1 V, Vd=0 V` using the CMD's nonlinear Poisson budget, then prepares
the Poisson/DG equilibrium. The CTest uses electrostatic potential and
referenced quasi-Fermi potentials with outer DG iteration. The verification
driver accepts `--dg-coupling monolithic_box` for simultaneous transport and
quantum Newton updates in the same quasi-Fermi coordinates.
The drain step bounds, extrapolation setting, and 21 requested
CurrentPlot voltages come from the CMD.

The fixture pins the DG drain current from Sentaurus T-2022.03-SP2
`IdVds_n13_des.plt`, as recorded in the historical comparison linked
above. It applies the existing 5% maximum relative current error over
0.1–2 V and 3% at 2 V. The report includes the zero-bias current and
four-terminal current sum. It also requires self-consistent DG convergence
at every accepted point and checks finite, nonzero nodal electron quantum
potentials in the 0 V and 2 V HDF5 states, including a response to the drain
bias. The field gate checks Vela's self-consistent DG state and its bias
response. The Sentaurus comparison uses terminal current.

## DG Id-Vg

Enable `VELA_ENABLE_TRANSPORTMODELS_DG_IDVG_CI=ON` to register
`reference_u2022_transportmodels_dg_idvg`. It uses the same grid, material
staging, oneMKL PARDISO backend and HDF5 state support as DG Id-Vd.
The standalone `pp7_des.cmd` retains the 300 K isothermal physics and the
original `Quasistationary` controls. It starts at `Vg=-1 V, Vd=0 V`, ramps
the drain to 1.1 V, then sweeps the gate to 2.2 V at 21 requested output points.
The two sweeps both consume the preceding stage's state.
The transport block uses potential and contact-basin quasi-Fermi
potentials with up to 80 Newton updates per frozen-quantum transport solve.
The DG inner Newton solve starts its line search at a full step. It still
backtracks on insufficient residual decrease and uses the existing inner
and outer convergence criteria.

The fixture pins `pp7_des.par`, the shared PolySi excerpt, and the native
T-2022.03-SP2 drain-current values from `IdVgs_n7_des.plt`. The manifest
records the original PLT hash and the archived comparison containing the
extracted values. The regression checks self-consistent electron quantum
potential convergence and the saved quantum fields across the three stages.
It compares published PLT currents in A/um with the existing Id-Vg
contract, including the 0.15 dex low-current limit, 10% on-state limit,
and the minimum drain-current-to-KCL-residual ratio of 10 at the three
deep-off points.
