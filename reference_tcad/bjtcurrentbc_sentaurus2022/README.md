# BJTcurrentBC single-base-current reference

This CTest uses the first base-current setting from the U-2022.12
Applications Library project `GettingStarted/sdevice/BJTcurrentBC`. It checks
conversion from voltage to current control at the base and the subsequent
collector sweep while base current remains at 1e-8 A/µm. The emitter stays at
0 V. The test compares the base voltage and the three terminal currents at 21
collector voltages from 0 to 1 V.

Generate the native inputs from the library's boundary TDR and mesh deck by
resolving its InputFile and running `snmesh --max_threads 1 bjt_msh`. Export
`bjt_msh.tdr` using `tdx -dd -M 0 -S 0 bjt_msh.tdr bjt`. Expand the library
SDevice template through the first base-current goal (Ib_ini = 1e-8), then run
one of the original collector sweeps using the immediately preceding state
and request 20 equal CurrentPlot intervals. Use the silicon defaults at
300 K. Keep the raw TDR, PLT and logs outside version control. The manifest
pins the actual native inputs and output and the neutral fixture files.

The Vela deck uses the same DF-ISE grid, dopant species and active models. Its
parameter excerpt pins the SRH and Auger overrides to the U-2022.12 Silicon
values. The Vela material defaults supply the remaining parameters. The test
uses the same U-2022.12 constants and oneMKL Release reference setup as the
other reference cases. The comparison includes base-current tracking and KCL
at every requested point. Collector and emitter current errors use a scale
floor equal to the base setpoint around zero crossings.

## Contact high-field mobility

The sidecar enables `solver.mobility.contact_electric_field_fallback`.
Sentaurus uses the electrostatic field magnitude in semiconductor elements
touching a contact for its default GradQuasiFermi high-field saturation.
Interior elements use the electron or hole quasi-Fermi gradient. See
the *Sentaurus Device User Guide*, chapter 15, equations 352–353.

The sidecar selects `element_box_vertex` for the Masetti high-field drive.
The solver applies saturation at each transport triangle's vertices before
averaging the mobility with nodal Box measures and matching edge coefficients.
The fixture uses Delaunay-transferred Box geometry, cell-material permittivity,
and signed transport volumes for Poisson charge. The contact-rule comparisons
below establish the earlier edge-based result.

On the 3,659-node, 7,064-triangle silicon mesh at 300 K, with
`Vb = 0.65 V` and `Vc = 0.5 V`, the native base current was
5.86233 nA/µm. Vela produced 3.89627 nA/µm before the contact rule and
5.87123 nA/µm with it. Explicit `ComputeGradQuasiFermiAtContacts = UseQuasiFermi`
in Sentaurus produced 3.87208 nA/µm. With high-field saturation off in both
programs, the native and Vela base currents were 5.91168 and 5.92036 nA/µm.

At fixed collector voltage, write `gB = ∂Ib/∂Vb` and `gC = ∂Ic/∂Vb`.
Keeping `Ib` constant requires the first-order bias correction
`ΔVb = -ΔIb/gB`. The resulting collector-current change is
`ΔIc - (gC/gB) ΔIb`, where the two current differences on the right are
measured at fixed voltage. This explains the larger collector-current
error in the original constant-base-current sweep.

On the same native `Vb = 0.5 V, Vc = 0 V` state, both Jacobians and RHS
vectors were compared on shared free rows after converting population
coordinates and normalizing each row by its own diagonal. For the 27
contact-neighbor nodes, the electron RHS difference L2 fell from 0.5510
to 0.0201 and the hole difference from 1.3819 to 0.3580. The interior
RHS differences changed by less than 1e-8 in these normalized units with
the original edge-projected interior drive.

## Linux Release verification

All four solve stages completed and all 21 requested collector points
passed using oneMKL PARDISO and the U-2022.12 physical constants profile.

| Maximum error | Before contact rule | With contact rule | CTest limit |
| --- | ---: | ---: | ---: |
| Base current tracking (A/µm) | 9.99e-13 | 9.68e-13 | 1e-12 |
| Base voltage (mV) | 11.22 | 0.0442 | 5 |
| Collector current (normalized) | 52.92% | 0.2112% | 5% |
| Emitter current (normalized) | 51.74% | 0.2041% | 5% |
| Terminal current imbalance (A/µm) | 8.38e-17 | 5.86e-17 | 1e-11 |

## Newton coordinates and step comparison

The case uses `potential_density` with `primary_variable_coordinates=unit`
throughout. The Newton variables are electrostatic potential and the two
carrier populations. The solver uses `ψ/VT`, `n/C0` and `p/C0`, and
exports volts and cm⁻³. Both voltage ramps use the CMD's normalized
`InitialStep=0.01`, `MinStep=1e-6`, `MaxStep=0.2`, `Increment=1.3`
and default `Decrement=2`, with Newton-work-dependent step growth.
The acceptance counts refer to steps past the initial point.

The current ramp in Vela solves a scalar current constraint at the requested
current targets. Sentaurus includes the contact equation in its coupled
Newton system. During the collector sweep, Vela adapts the outer voltage
using the last device solve's iterations. The CSV separately reports the
number of device solves and total Newton updates within each successful
current root, and VTU metadata reports the solved base voltage at each point.

On the output-instrumented comparison, both programs used six Poisson
updates and one coupled equilibrium update. The base voltage ramp accepted
14 steps per program and used 37 native and 38 Vela Newton updates. The
collector ramp accepted 28 native steps and 44 Vela steps. Its native
accepted steps used 98 Newton updates. The last Vela device solve at each
accepted step used 135 updates in total, and the successful current roots
used 91 device solves and 283 updates.

A separate current-ramp replay reached the 15 normalized current targets
from the native accepted trajectory. It took 61 Vela device solves and
119 Newton updates inside the current roots. The native coupled solves
used 49 updates at these targets.

The field comparison matched all 3,659 nodes at each of the 21 collector
output points. The maximum potential difference was 9.92 mV, with a
maximum pointwise RMS of 0.483 mV across the 21 bias points. The largest
per-bias 95th-percentile concentration errors were 22.01% for electrons
and 2.86% for holes, using `max(|Sentaurus|, 1e7 cm^-3)` in the denominator.
The pointwise maxima were 194.65% and 31.38%, respectively.

Configure with `VELA_ENABLE_BJTCURRENTBC_CI=ON` to register the test.

```sh
ctest --test-dir <build> -R '^reference_u2022_bjtcurrentbc$' --no-tests=error --output-on-failure
```
