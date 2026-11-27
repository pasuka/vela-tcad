# Reference TCAD regression tests

`poisson_arithmetic_probe` runs the built CLI against a small semiconductor
fixture, checks the read-only arithmetic export against production binary128
residuals, preserves inputs and existing output files, and rejects unsupported
Fermi statistics. The numerical material/Poisson tests also check sub-ULP charge
and dielectric responses against the physical Jacobian. These tests do not
qualify a full split-state IALMob/Auger/avalanche Newton path.

`python -B tests/regression/test_slot_newton_stabilization.py` checks that the
two independent SLOT stabilization controls preserve physics and acceptance
settings, reject invalid enabled scales, and do not overwrite existing evidence.
Numerical four-arm evidence is recorded separately in the validation reports.

The SLOT IdVd driver tests (`test_run_slot_ldmos_idvd.py`) cover the optional
`--gate-step-policy adaptive --gate-only` route, including the Fermi Poly
predictor advisory monitor. `GATE_BIAS_WARNING_THRESHOLD_V=0.1` writes a warning
for every larger bias change; it does not reject or clip the prediction. The
adaptive policy retains its separate 0.25 V maximum. Neither value is a device
acceptance tolerance: prediction must still pass Poisson preparation and full
Newton qualification. Tests cover repeated warnings, continuation from a
qualified nonzero gate bias, and halving the actual failed span. An optional
explicit cap in the generic policy is also tested, but is not applied by the
Fermi gate route.
Mocked continuation to 5 V establishes control flow only, not numerical
qualification. Existing frozen campaign directories are not updated by changes
to the source driver.

This directory contains Python regression tests for checked-in cross-TCAD
fixtures, import/conversion tools, comparison scripts, and validation report
contracts.  Device-level reusable inputs live under `reference_tcad/`; the
repository no longer ships a separate set of uncalibrated engineering example
decks.

The Python regression environment needs NumPy, Pillow, h5py, and jsonschema. On Ubuntu,
install `python3-numpy python3-pil python3-h5py python3-jsonschema` and use `/usr/bin/python3`
for CMake/CTest so that it sees the apt-installed modules. The C/C++ HDF5
development library alone does not provide the Python `h5py` module.

SDE modeling regressions additionally need the external Gmsh 4.15.2 Python
wrapper, shared library and CLI. The jobs in `.github/workflows/ci.yml` prepare them with
`scripts/prepare_sde_external_gmsh.py` using the versions and SHA-256 hashes in
`configs/sde_external_runtimes.json`; this exports `VELA_GMSH_PYTHON`,
`VELA_GMSH_LIBRARY` and `VELA_GMSH_EXECUTABLE` to later CI steps. On Ubuntu,
install `libglu1-mesa python3-venv` for that runtime and `patchelf` for the
Linux package relocation regression and offline package assembly. The
apt-installed `python3-jsonschema` supplies its transitive Python dependencies;
the isolated offline package instead uses the six pinned wheels in
`configs/sde_offline_inputs.json`. No additional C++ link dependency is required
for this modeling path.

For local Windows regression builds, install the UCRT64 `python-numpy`, `python-pillow`, and
`python-h5py` and `python-jsonschema` packages and pass
`-DPython3_EXECUTABLE="$(cygpath -m /ucrt64/bin/python.exe)"` to CMake from the
MSYS2 shell. Check imports with `/ucrt64/bin/python.exe` as well: an unpinned
CMake search can select the hosted runner's separate Python installation even
when the shell's `python` is UCRT64.

On Windows, a `h5py.defs` DLL "procedure not found" error can also indicate a
feature mismatch within the same HDF5 DLL name. The
[ROS3 import repair record](../../docs/validation/h5py_ros3_dll_repair_2026-10-25.md)
documents an h5py package importing ROS3 symbols from an HDF5 build without
ROS3. Rebuild/install the binding against the actual HDF5 feature set; preserve
the import failure as a failure until repaired. Do not turn it into a test skip.

Sentaurus frontend tests require full JSON Schema validation, including nested
properties. Missing `jsonschema` is an environment failure, not a reason to
silently fall back to top-level structural checks.

LDMOS linked-input manifests retain byte-exact SHA-256 checks. The repository's
`.gitattributes` specifies the originally qualified LF or CRLF form for each
checked-in dependency; do not replace hashes or normalize bytes in the digest
function to make a platform-specific mismatch pass.

The `ialmob_element_avalanche_probe_protocol` test checks the diagnostic
`ialmob_element_avalanche_probe` SI contract, explicit source shares, output
derivative columns and rejection of unsupported temperatures/RefDens. Full local
chain-rule differences (including live mobility and both carrier statistics) are
in `test_impact_ionization` under `[element-density]`. These probe checks do not
qualify coupled SLOT avalanche or BVDS continuation.

Run the main reference-tool checks from the repository root:

```bash
python -m unittest tests.regression.test_reference_tcad_tools
ctest --test-dir build --output-on-failure -R reference_tcad_regression
```

Many reference cases have additional focused test modules in this directory.
Each case README or machine-readable `*_reference.json` inventory identifies
its authoritative scripts, reports, and acceptance boundary.

## PN2D Release regression

The `CI` workflow includes PN2D in the full Release CTest run on every push and
pull request. To run only PN2D manually, open **CI** under **Actions**, choose
**Run workflow**, select the branch and set `suite=pn2d`. Choose `all` to include
the full Release suite and SDE offline package qualification.
The default branch must expose the workflow before GitHub offers manual runs.
With GitHub CLI, use `gh workflow run ci.yml --ref <branch> -f suite=pn2d`.

PN2D uses the shared Release build, `sentaurus_u2022_12` constants and the checked-in
[reference archive](../../reference_tcad/pn2d_sentaurus2018/ci/README.md).
It runs equilibrium, the complete IV and compensated BV sweeps, and both
compensated sweeps for all eight M3 variants through `sd-sim`. It checks the
requested biases, current and conservation gates, and saved spatial fields.
The staged CI sidecar sets `stage_vtk_bias_points` for spatial comparison and
the final device plot. A VTU file is expected for each selected bias. The
solver still runs the complete CMD sweep and checks all required current
points. See the [SDevice output contract](../../docs/sentaurus_sdevice_workflow.md).
The compensated BV and all sixteen M3 sweep CTests select UMFPACK when CMake
finds it. Builds without UMFPACK run those same cases and acceptance checks
with Eigen SparseLU.
The SDE test generates and reimports the PN2D mesh twice with pinned Gmsh.
The job uploads JSON reports and logs as `pn2d-release-regression`.

For a local Linux run with the build dependencies and Gmsh runtime above:

```bash
cmake -S . -B build-pn2d -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DPython3_EXECUTABLE=/usr/bin/python3 \
  -DVELA_ENABLE_PN2D_CI=ON -DVELA_PHYSICAL_CONSTANTS_PROFILE=sentaurus_u2022_12
cmake --build build-pn2d --target sd-sim vela_example_runner --parallel 2
ctest --test-dir build-pn2d --output-on-failure --no-tests=error -L '^pn2d-ci$'
```

`sdevice_e2e_gates` runs five Schottky and BVmethods deck checks from
[pinned TDX GRD/DAT fixtures](../../reference_tcad/sdevice_check_only/manifest.json).
It also checks rejection of the original TransportModels deck with thermal
boundaries. The isothermal DD Id-Vd adaptation has a separate numerical test.
The test uses `vela device --check` and requires no external reference directory.

`simplemos_evidence_archive` checks the selected historical evidence, archive
manifest and safe restoration with temporary fixtures. It needs no VM access;
full external archive verification is a separate explicit command documented in
[the recovery guide](../../reference_tcad/simplemos_sentaurus2022/local_evidence_20260930/README.md).

Generated TDR, VTK, accepted-state, and log files remain under ignored
`build*/reference_tcad/` directories.  Tests should use checked-in neutral
inputs or explicitly generated temporary fixtures and must not present
synthetic smoke data as commercial-tool calibration evidence.

SimpleMOS historical M9/M10/M12 and M30-M46 checks validate original report
and artifact bytes against archived source provenance, not the evolving live
source tree. See `tests/fixtures/simplemos_historical_sources/README.md` for
the 132 recovered source revisions and the one explicitly unrecovered source
hash shared by M33/M35/M37. The archive-integrity test asserts this gap; it does
not imply complete historical reproducibility. Current solver qualification
uses numerical/property tests and a separately sealed engineering campaign.

The opt-in SLOT neighboring-row recovery properties are covered by
`test_slot_local_row_repair.py` in the `slot_ldmos_loadline_continuation` CTest:
disabled mode performs no work, patch selection preserves all violating rows,
connected weak neighbors remain in the patch, backtracking checks spillover
across the full active row set, and row descent cannot override the original
global screen. Unscreened trials are never selected. The separate
`--predict-boundary` option ranks same-carrier neighboring rows by finite-difference
response, expands in the first cycle within the existing eight-DOF cap, and reuses
already evaluated columns; it does not change any acceptance threshold.
Synthetic callback tests verify
proposal logic only; real-device qualification still requires the unchanged
Newton solve, full row audit and bound evidence.

The retired custom NMOS, coarse7x3, Minimal6, and skewed-Tri3 reference cases
are no longer regression inputs. Import/contact and mesh invariants remain
covered using temporary data and retained fixtures. Shared avalanche helper
tests are `test_sentaurus_avalanche_replay.py` and
`test_sentaurus_avalanche_controls.py`; the standalone contact override test
is `test_sentaurus_contact_overrides.py`.

The `slot_ldmos_noavalanche_continuation` CTest entry includes guarded Poly
restart helper tests. For the SLOT continuation script, the opt-in
`--restore-poly-equilibrium` flag requires carrier corrector passes and a local
continuity audit. It prepares isolated source-free PolySilicon candidates only;
full coupled convergence, final all-active row checks, and update-history
qualification remain mandatory. Unsupported profiles are rejected rather than
silently projected. See the [2026-09-22 SLOT evidence](../../docs/validation/slot_ldmos_high_voltage_resume_2026-09-22.md).

`slot_ldmos_idl_configuration` checks the independent IDL Bennett-Wilson
material bundle, a cold zero-bias start, and rejection of incomplete continuity
proof. Run directly with `python -m unittest tests.regression.test_slot_idl_configuration`.
The numerical path is deliberately separate: `scripts/run_slot_ldmos_idl.py`
requires the archived SLOT mesh/native geometry resources and a specified runner.
It starts without a saved state, establishes Vd=0.1 V, then sweeps Vg to 5 V.
Density/poly restart operations provide candidates only; acceptance requires
the original complete Newton, active-row and global-closure gates. This IDL
material bundle must not silently replace the historical BVDS materials.

The separate IdVd package is checked with
`python -m unittest tests.regression.test_slot_idvd_configuration`.
`scripts/run_slot_ldmos_idvd.py` starts cold, establishes Vg=5 V at Vd=0,
then follows the supplied native drain samples with at most 1 V intermediate
steps. It freezes resources and the runner, bounds child and total runtimes,
halves a failed interval at most three levels, and accepts only original
full-device gates. It deliberately does not reuse the Boltzmann IDL poly reset.
Its `goal_reached` means completion of the qualified voltage path; current
agreement with Sentaurus is separately reported at identical drain biases.
`test_run_slot_ldmos_idvd.py` checks that a successful process with a rejected
device state cannot become a retry parent; it uses mocked rows and is not
numerical qualification.

The opt-in `--drain-initialization secant` uses the two most recent qualified
drain states, with a target/history step ratio in (0,3]. The default remains
`warm`. A failed secant solve gets one same-target warm attempt before the
existing bounded bisection; an out-of-range ratio uses warm initialization.
Rejected states never enter either parent history. The Fermi split-QF/log-density
predictor is an unqualified guess and does not change the constitutive model.
`--drain-import-summary` checks the last two saved states of a completed batch
against full-gate IV records, hashes, physical inputs and Vg=5 V contact metadata.
Imported records are separate from new records; this check is not an independent
numerical re-solve. `--drain-stop` limits the requested segment to an exact native
sample (default 60 V). For imported/short runs, `goal_reached` refers only to that
requested segment, as recorded by `requested_drain_segment_qualified`.
Run the corresponding control, predictor and configuration tests with:

```text
python -m unittest tests.regression.test_run_slot_ldmos_idvd tests.regression.test_slot_idvd_secant_predictor tests.regression.test_slot_idvd_configuration
```

The explicit `--gate-initialization fermi_poly_poisson` option adds an isolated
Poly Fermi bias predictor and a `preserve_initial` Poisson stage before each
gate-point full Newton solve. The default remains `none`. Predictor/Poisson
states never enter the accepted records or become retry parents, and the drain
stage retains full transport Newton without this equilibrium initialization.
All substage time and failures count toward the same campaign budget. The
driver regression covers the handoff, rejection isolation and gate-only scope.

`python -m unittest tests.regression.test_slot_fermi_poly_predictor` checks the
explicit Fermi Poly rigid-bias initial guess: preserve occupation arguments,
densities and unselected nodes, and reject unsupported packed coordinates.
`scripts/slot_fermi_poly_predictor.py` does not qualify a state or change the
production IdVd driver. Its output requires Poisson/coupled requalification;
it does not apply a Boltzmann density reset to a Fermi configuration.

The `slot_ldmos_loadline_continuation` test covers the SLOT-specific independent
device load-line driver, `scripts/run_slot_ldmos_loadline_continuation.py`.
It freezes R=1e12 ohm*um, a 1e-6 V load-line tolerance, and the original SLOT
all-active-row/global-closure gates. `--initial-inner-bias` is the device voltage;
`--outer-targets` are supply voltages. A seed is requalified before use. Each new
device bias needs actual coupled Newton updates, and each external acceptance
requires byte-identical circuit confirmation. A changed confirmation state stops
incomplete; it cannot inherit the previous proof.

Choose a short, fresh ignored `--output-dir`. `--resume` requires identical input
files, solver/helper code, runner, solver/thread environment, request, and budgets.
It validates saved evidence and does not repeat completed solves. Numerical
failure permits bounded inner/outer step reduction; missing runtime files,
corrupt proof, or resource/environment changes stop the campaign. This entry is
qualified only for the documented SLOT profile, not a general circuit solver or
proof of original AvalDensGradQF parity. See the
[2026-09-23 review](../../docs/validation/slot_ldmos_loadline_review_2026-09-23.md).

The Catch2 `test_split_dd_runtime` target includes `[split_sg_diagnostics]`
properties: equilibrium with unequal carrier-specific population prefactors,
sub-ULP driving forces, edge-orientation reversal, common-ni fallback and explicit
unavailability for edges without transport support. Diagnostics must leave the
production residual, current and low coordinates unchanged. Its 200-digit SG
reference holds production mobility fixed and does not qualify a device state.

### BVMethods full breakdown CI

Enable `VELA_ENABLE_BVMETHODS_CI=ON` in a Release build with
`VELA_PHYSICAL_CONSTANTS_PROFILE=sentaurus_u2022_12`. Build `sd-sim` and
`vela_example_runner`, set `VELA_MKL_RT_LIBRARY` to the installed oneMKL runtime
and expose its dependencies through the platform library search path. Run
`ctest --test-dir <build> --output-on-failure --no-tests=error -L "^bvmethods-ci$"`.

The CTest fixture stages the SHA-256-pinned DF-ISE pair in
[`reference_tcad/sdevice_check_only`](../../reference_tcad/sdevice_check_only/manifest.json).
For a fresh export with the qualified TDX and original TDR, configure
`-DVELA_BVMETHODS_TDX=<tdx-executable>` and
`-DVELA_BVMETHODS_TDR=<n1_msh.tdr>`. The fixture calls
`tdx -dd -M 0 -S 0 n1_msh.tdr n1_msh` and checks the output hashes against the
same manifest. The GitHub runner uses the pinned files directly.

The test runs the qualified U-2022.12 BVMethods CMD and sidecar through
`sd-sim` with single-thread PARDISO. It checks the stage chain, the original
drain current stop, the interpolated BV, the native 0-100 V current curve and
the accepted-point numerical checks. The native curve is bounded at 100 V. Its interpolated current gate is 2%
with a 1e-12 A/um denominator floor. BV uses the first 1e-4 A/um
crossing and a 0.5% tolerance against the qualified 6.673833 V anchor.
Reports are written to `<build>/bvmethods-reports/`.

The shared fixture also supplies `bvmethods_ci_ABA_poisson`,
`bvmethods_ci_ABA_coupled`, and `bvmethods_ci_voltage2current`. The versioned
[method references](../../reference_tcad/bvmethods_sentaurus2018/u2022_12_20261117_methods/reference.json)
use native U-2022.12 runs on the same mesh and material inputs, with BTBT off.
The two ABA regressions check the physical integral stop and its crossing
bracket. The switched run checks the 6 V handoff, current tracking, KCL and
current-interpolated BV. Reports use separate subdirectories for the three
methods.

## SDevice numerical reference regressions

Enable `VELA_ENABLE_REFERENCE_U2022_CI` to register `reference_u2022_gate_contacts`,
`reference_u2022_schottky_iv`, `reference_u2022_bjt_m0`, and `reference_u2022_bjt_m1`.
Use `ctest --test-dir <build> -N -L '^reference-u2022$'` to check registration
for the configured source revision. These tests use the built `sd-sim` with
the U-2022.12 constants profile. CI enables the group in its full suite and
exposes `suite=reference-u2022` for a manual run.

Gate checks zero barrier, +0.55 V and -0.30 V barriers, and both zero and
nonzero work-function offsets against the frozen analytical flat-band references.
It checks 41 surface potential snapshots per CMD and brackets the flat-band
crossing. The flat-band voltages, signed barrier shifts and zero work-function
shift must each agree with their reference values within 1 mV. The
[Gate manifest](../../reference_tcad/gate_reference/ci_u2022_12/manifest.json)
pins a GRD/DAT mesh generated from `gate_sde.cmd` using Gmsh 4.15.2.

The [Schottky manifest](../../reference_tcad/schottky_charon_sentaurus2018/ci_u2022_12/manifest.json)
and [BJT manifest](../../reference_tcad/genius_bjt_sentaurus2022/ci_u2022_12/manifest.json)
pin the tested CMDs, parameter overrides, numerical sidecars and normalized
reference currents. Schottky reuses the existing TDX GRD/DAT export. The BJT
GRD/DAT files serialize the validated neutral mesh and species concentrations.
Checksums use LF line endings. The tests run without Sentaurus, TDR or Gmsh.

Schottky checks the complete 0-to-1 V forward sweep against 24 native
positive-bias points (0.5 dex), monotonicity, and terminal current balance.
BJT M1 initializes at 0 V, ramps the base to 0.70 V, and compares both terminal
balance and all 31 collector output points through 3.0 V. In the 0.5-to-3.0 V
window, each terminal current and the current gain must agree to 0.05 dex.
Its numerical regression uses EFT for the linear defect checks and records
the requested collector points in CSV. `sd-sim` publishes only the final
collector field, with intermediate runner artifacts in the report directory.

BJT M0 uses the same GRD/DAT mesh and retains Boltzmann statistics and SRH.
Its CMD explicitly selects Bennett-Wilson, with coefficients from the archived
T-2022.03-SP2 PAR, to reproduce the native M0 default.
It runs through `sd-sim` in Unit potential-density coordinates and requires
the 0.70 V base stage, all 31 collector points, finite terminal currents,
active-region current signs, and the existing KCL limits. It reports current
and gain differences against the archived T-2022.03-SP2 M0 curve for
characterization. The M1 0.05 dex gate applies to its own U-2022.12 reference.

Enable `VELA_ENABLE_TRANSPORTMODELS_IDVD_CI` to register the TransportModels
qualification test while the full-sweep solver is being validated. It compares
the isothermal DD Id-Vd sweep at Vg = 1 V against 21 archived T-2022.03-SP2
output points. See the
[TransportModels fixture](../../reference_tcad/transportmodels_sentaurus2022/ci_u2022_12/README.md)
for the input provenance, step controls and current acceptance criteria.

The DD Id-Vg test `reference_u2022_transportmodels_idvg` is enabled by default
for fresh `sentaurus_u2022_12` builds with oneMKL PARDISO. Its 21 comparison
points use unit-scaled electrostatic potential and carrier concentrations,
with compensated state arithmetic on the three deep-off points.

### DG reference tests

The CMake options below are off by default for local configurations. The
full Release configuration in the `CI` workflow sets them as shown.

| Test | CMake option | Full Release CI | DG coupling used by CTest |
| --- | --- | --- | --- |
| `reference_u2022_transportmodels_dg_idvd` | `VELA_ENABLE_TRANSPORTMODELS_DG_IDVD_CI` | On | `monolithic_box` |
| `reference_u2022_transportmodels_dg_idvg` | `VELA_ENABLE_TRANSPORTMODELS_DG_IDVG_CI` | Off | `outer` when explicitly enabled |
| `reference_u2022_singledevice_idvg_linear` | `VELA_ENABLE_SINGLEDEVICE_IDVG_CI` | On | `monolithic_box` |

The two monolithic selections are set by CTest. See the
[TransportModels fixture](../../reference_tcad/transportmodels_sentaurus2022/ci_u2022_12/README.md)
and the [SingleDevice fixture](../../reference_tcad/singledevice_sentaurus2018/ci_u2022_12/README.md)
for the bias ranges and acceptance gates.

### Running the U-2022.12 group locally

BJT and TransportModels use single-threaded oneMKL PARDISO. Set `VELA_MKL_RT_LIBRARY` to the
absolute path of the installed MKL runtime. Schottky uses UMFPACK when
available and Eigen SparseLU otherwise, as does Gate. Configure and run the opt-in group
with a Release build:

```bash
cmake -S . -B build-u2022 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DVELA_PHYSICAL_CONSTANTS_PROFILE=sentaurus_u2022_12 \
  -DVELA_ENABLE_MKL_PARDISO=ON -DVELA_ENABLE_REFERENCE_U2022_CI=ON
cmake --build build-u2022 --target sd-sim vela_example_runner
ctest --test-dir build-u2022 --output-on-failure --no-tests=error -L '^reference-u2022$'
```

Each run writes its configs, curves, logs and `result.json` to a fresh directory
under the build's `reference-u2022-reports/` directory.

### SDS numerical case labels

The PN2D, BVMethods and U-2022.12 registrations covered by the SDS
classification receive one device label when their CMake options are enabled.
The groups are PN2D (`sds-device-pn2d`), MOS / Gate (`sds-device-mos`), Schottky
(`sds-device-schottky`), and BJT (`sds-device-bjt`).

Preview the selected tests in an already configured build. CTest adds the
fixture setup required by the selected numerical tests.

```powershell
ctest --test-dir build-release/sentaurus2022 -N -L "^sds-device-pn2d$"
ctest --test-dir build-release/sentaurus2022 -N -L "^sds-device-bjt$"
```

The [desktop catalog](../../docs/vela_workbench.md#packaged-cases-and-ctest)
lists the packaged inputs. The BJT label also selects the source-only
BJT M0 regression. Select the two SLOT LDMOS tests with `slot-ldmos`.

## SLOT LDMOS device regressions

Configure with `VELA_ENABLE_SLOT_LDMOS_CI=ON` to register
`slot_ldmos_bvds_sdevice`, `slot_ldmos_idvg_sdevice` and
`slot_ldmos_idvd_sdevice`. The option defaults to `OFF`. Set
`VELA_SLOT_LDMOS_BVDS_CMD`, `VELA_SLOT_LDMOS_IDVG_CMD` and
`VELA_SLOT_LDMOS_IDVD_CMD` to preprocessed original CMD files in
self-contained input directories. Include the TDR or GRD/DAT mesh and
doping data, PAR files, and any required SLOT sidecar or restart state.
Use paths relative to the CMD for inputs and outputs. The device inputs
are supplied separately from the repository. If the IdVd package uses an
explicit material defaults library, pass its path as
`VELA_SLOT_LDMOS_IDVD_DEFAULT_PAR`.

IDSAT has a separate switch, `VELA_ENABLE_SLOT_LDMOS_IDSAT_CI`, which also
defaults to `OFF`. Enable it to register `slot_ldmos_idsat_sdevice` and pass
the prepared original CMD as `VELA_SLOT_LDMOS_IDSAT_CMD`. If the package
uses a separate defaults library, set `VELA_SLOT_LDMOS_IDSAT_DEFAULT_PAR`.
The three other SLOT CMD packages are optional when enabling IDSAT alone.

```bash
cmake -S . -B build-release \
  -DVELA_ENABLE_SLOT_LDMOS_IDSAT_CI=ON \
  -DVELA_SLOT_LDMOS_IDSAT_CMD=/path/to/idsat/IDSAT_des.cmd \
  -DVELA_SLOT_LDMOS_IDSAT_DEFAULT_PAR=/path/to/idsat/defaults.par
cmake --build build-release --target sd-sim
ctest --test-dir build-release -R '^slot_ldmos_idsat_sdevice$' --output-on-failure
```

Build `sd-sim` and its runtime targets, then run `ctest -L slot-ldmos --output-on-failure`.
The runner uses the build's `sd-sim` by default. `VELA_SLOT_LDMOS_SD_SIM`
can explicitly select a separately qualified SLOT-capable runtime with its
companion tools and scripts. The tests require successful `sd-sim`
execution and verify potential and density unknowns with UnitScaling in
every emitted stage. With `Math { ExtendedPrecision }`, each stage must
report end-to-end float80 and produce a verified accepted-state archive.
The tests check the CMD voltage endpoint or the BVDS current stop, as well
as the current curve, final VTU and run log. IdVd must first complete the
gate sweep at zero drain bias before sweeping the drain. IDSAT checks both
voltage endpoints and verifies that the drain endpoint at zero gate bias
is saved and reloaded before sweeping the gate at fixed drain bias. Each test keeps
its run directory and logs under `slot-ldmos-reports/` in the build tree.
