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
wrapper, shared library and CLI. Both jobs in `.github/workflows/ci.yml` prepare them with
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
