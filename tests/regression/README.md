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

On Windows CI, install the UCRT64 `python-numpy`, `python-pillow`, and
`python-h5py` and `python-jsonschema` packages and pass
`-DPython3_EXECUTABLE="$(cygpath -m /ucrt64/bin/python.exe)"` to CMake from the
MSYS2 shell. Check imports with `/ucrt64/bin/python.exe` as well: an unpinned
CMake search can select the hosted runner's separate Python installation even
when the shell's `python` is UCRT64.

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
