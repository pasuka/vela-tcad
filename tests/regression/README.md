# Reference TCAD regression tests

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
The staged CI sidecar selects the VTU biases used for spatial comparison and
the final device plot. All CMD voltage targets remain in the sweep.
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

The retired custom NMOS, coarse7x3, Minimal6, and skewed-Tri3 reference cases
are no longer regression inputs. Import/contact and mesh invariants remain
covered using temporary data and retained fixtures. Shared avalanche helper
tests are `test_sentaurus_avalanche_replay.py` and
`test_sentaurus_avalanche_controls.py`; the standalone contact override test
is `test_sentaurus_contact_overrides.py`.

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
`reference_u2022_schottky_iv`, and `reference_u2022_bjt_m1`. They run their
CMDs through the built `sd-sim` with the U-2022.12 constants profile. CI enables
them in the full suite and exposes `suite=reference-u2022` for a manual run.

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
characterization. The U-2022.12 M1 parity gate retains its own reference.

Enable `VELA_ENABLE_TRANSPORTMODELS_IDVD_CI` to register the TransportModels
qualification test while the full-sweep solver is being validated. It compares
the isothermal DD Id-Vd sweep at Vg = 1 V against 21 archived T-2022.03-SP2
output points. See the
[TransportModels fixture](../../reference_tcad/transportmodels_sentaurus2022/ci_u2022_12/README.md)
for the input provenance, step controls and current acceptance criteria.

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

Numerical regressions that explicitly launch `sd-sim` receive one SDS device
label when their existing CMake options are enabled. The groups are PN2D
(`sds-device-pn2d`), MOS / Gate (`sds-device-mos`), Schottky
(`sds-device-schottky`), and BJT (`sds-device-bjt`).

Preview the selected tests in an already configured build. CTest adds the
fixture setup required by the selected numerical tests.

```powershell
ctest --test-dir build-release/sentaurus2022 -N -L "^sds-device-pn2d$"
ctest --test-dir build-release/sentaurus2022 -N -L "^sds-device-bjt$"
```

## SLOT LDMOS device regressions

Configure with `VELA_ENABLE_SLOT_LDMOS_CI=ON` to register
`slot_ldmos_bvds_sdevice` and `slot_ldmos_idvg_sdevice`. Set
`VELA_SLOT_LDMOS_BVDS_CMD` and `VELA_SLOT_LDMOS_IDVG_CMD` to preprocessed
CMD files in self-contained input directories. Include the input mesh,
doping, parameter files, and any SLOT sidecar and restart state required by
the CMD. Use paths relative to the CMD for inputs and outputs. The device
inputs are supplied separately from the repository.

Build `sd-sim` and its runtime targets, then run `ctest -L slot-ldmos --output-on-failure`.
The runner uses the build's `sd-sim` by default. `VELA_SLOT_LDMOS_SD_SIM`
can explicitly select a separately qualified SLOT-capable runtime with its
companion tools and scripts. The tests require successful `sd-sim`
execution, verify potential and density unknowns with unit scaling in the
emitted stages, and check the completed voltage sweep or BVDS current stop.
Each test keeps its run directory and logs under `slot-ldmos-reports/` in
the build tree.
