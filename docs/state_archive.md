# Restart state storage

Production restart files use `.h5` and the `vela.state/2` HDF5 schema. HDF5 and
HighFive are required by the standard build. The optional Sentaurus TDR importer
is a separate feature; a Vela archive is not a TDR file.

## Disk contract

The root attributes are `schema` (string), `node_count` (scalar uint64) and
`metadata_json` (UTF-8 JSON). Each `/fields/<name>` dataset is a one-dimensional,
uncompressed IEEE binary64 little-endian array of `node_count` elements. Its
`unit` attribute is mandatory. No text round trip is needed to read numerical
values into the solver or predictor.

| Fields | Unit | Requirement |
|---|---|---|
| `psi`, `phin`, `phip` | V | All states |
| `electrons_m3`, `holes_m3` | m^-3 | DD requires both; electrothermal allows both or neither |
| `electron_qf_reference_V`, `hole_qf_reference_V`, `electron_qf_increment_V`, `hole_qf_increment_V` | V | All four together when split coordinates are present |
| `electron_quantum_potential_V` | V | When available |
| `electron_quantum_potential_like_V` | V | Requires the quantum-potential field |
| `temperature_K` | K | Required in electrothermal mode; absent in DD mode |
| `population_electrons_internal`, `population_holes_internal` | internal_concentration | Optional complete potential-density DD coordinate pair, with explicit scale metadata |
| `packed_psi`, `packed_electron_qf_increment`, `packed_hole_qf_increment` | 1 | Optional complete normalized DD coordinate bundle; requires all four QF reference/increment fields |
| `packed_psi_low`, `packed_electron_qf_increment_low`, `packed_hole_qf_increment_low` | 1 | Optional complete low-component bundle; requires the high-component bundle |

The four-equation restart unknowns are potential, the two QF coordinates and
temperature. Density is derived by the physics implementation; the archive does
not invent zero densities to fill a DD-shaped record. Physical QF values and
split reference/increment values are checked for consistency. Keeping the split
arrays separately preserves increments smaller than the ULP of the reference.

The optional packed DD bundle records `packed_potential_scale_V` in metadata.
When low components are present, metadata also contains
`split_state_schema: vela.split-dd-state.v1` and `split_mesh_fingerprint`.
These are the explicit split-DD solver coordinates, distinct from the ordinary
reference/increment QF representation. Readers reject partial bundles, invalid
scales, and disagreement with the physical fields. The DD adapter restores all
three high/low blocks without a physical-potential round trip. Low tails,
including subnormals, retain their binary64 bits. The solver additionally checks
the split mesh fingerprint when loading its runtime; the archive still requires
the independent full `mesh_sha256` check.

The explicit classical IALMob split profile consumes this same complete bundle
in Poisson, populations, SG, recombination and local field/source evaluation.
Read-only residual packing and terminal evaluation preserve the checkpoint's QF
reference frame; they must not select a new terminal-nearest frame from the
final physical fields or replace it with geometric contact basins. Continuation
rebasing remains explicit and applies the new boundary conditions separately.
Public physical arrays
remain rounded projections, not substitutes for the packed high/low fields.
Supported physics and numerical limitations are listed under `split_dd_state`
in [the configuration schema](config_schema.md). See the
[BVDS production validation](validation/slot_ldmos_bvds_C_production_low_state_2026-10-16.md)
for full qualification, reload and explicit low-zeroing controls.
For `state_variables: potential_density`, the optional internal carrier pair
requires `population_coordinate_schema: vela.population-coordinates.v1` and a
positive finite `population_concentration_m3_per_internal` conversion factor.
Each internal value must be finite and nonnegative, and its binary64 product
with that factor must equal the corresponding SI density. A complete pair and
both metadata fields are required together. This bundle rejects split/packed
QF state. Physical `phin`/`phip` remain derived outputs; they do not replace n/p.
C++ and Python enforce the same storage contract.

The DD writer records the original carrier coordinates without a unit round
trip, including subnormals and inactive zeros. Same-scale restore reads these
coordinates directly. Different-scale restore converts SI output fields; older
population archives without this bundle remain readable through the SI path.
`DDSolution::populationRestoreMode` reports `exact_internal`, `converted_units`
or `legacy_si` respectively. The marker describes carrier restoration, not a
universal bitwise guarantee for potential renormalization or changed models.
The metadata `authoritative_fields` lists ψ and the internal carrier pair;
validation uses the explicit schema, field names and scale rather than trusting
that informational list. No cached coordinate copy is retained in DDSolution.

The migration-only CSV importer extracts repeated scale/schema/fingerprint
columns into metadata, checks their uniformity, and verifies every decoded
numeric field. Production continues to accept only HDF5. Python frame translation
rejects packed coordinates until an explicit transformation of all coordinates
has been supplied; it cannot silently translate only the physical fields.

Required metadata: `mode` (`dd` or `electrothermal`), `mesh_sha256` and finite
`potential_origin_V`. The canonical `vela.mesh/1` identity includes node order,
coordinates, cells, region/material names, contacts and length units. Readers
must supply an independently obtained expected identity and node count.
`state_archive.inspect()` is a reporting utility, not independent mesh validation.

Production writers also record source-configuration and physical-input file
hashes. DD output records `bias_V`, `bias_contact`, `contact_biases_V` and
`state_role`. DD bias values use the solver reference frame; adding
`potential_origin_V` recovers physical terminal voltage. Explicit frame
translation changes both fields and corresponding bias metadata. Predicted
states record their parent hashes and target bias. Electrothermal records carry
the sweep bias and source/boundary identity supplied by the sweep or point entry.

Source hashes are prepared before state writes. Sequential DD preparation reuses
them only after comparing input bytes, including additional mobility file bytes.
This does not reuse carrier values across temperature or change symbolic-analysis
invalidation. Configuration hashes include solver/model settings; they are
provenance, not a promise of cross-model qualification.

## Loading and checkpoint recovery

### Explicit SLOT reference-frame continuation

`scripts/run_slot_ldmos_reference_frame.py` prepares a **new** bundle with an
explicit `--target-origin-V`. Its `--initial-drain-voltage-V`, optional
`--previous-drain-voltage-V`, and `--targets` are absolute physical drain
voltages. These differ from the underlying direct driver's solver-frame
`--initial-bias`/`--targets`. Physical voltage equals stored voltage plus origin;
Vds is drain voltage minus source voltage in either frame.

The entry requires explicit Newton, enforced carrier-row and global-closure
gates, direct device bias, the SLOT Si/PolySilicon and SiO2/Oxide/Nitride material
classification, and ohmic contacts. External circuits, explicit boundary
objects and other contact/material profiles are rejected. This is not automatic
origin selection or a general device translation API.

The input config declares the current seed's source origin. A previous seed in
another frame requires `--previous-source-origin-V`; otherwise it must match
the input config's origin. Complete archive `bias_V`, `bias_contact` and
`contact_biases_V` provenance must match the declared physical voltages. Legacy
migration archives lacking these fields are rejected, rather than inferring
their bias from filenames. Mesh identity and origin are independently checked.

The helper preserves the compensated split-QF remainder and density arrays,
shifts potential/contact coordinates, and projects inactive contact dummy QF
constraints to the new terminal coordinates (62 SLOT nodes on the current
mesh). Interior inactive QF remains zero. Density arrays are input storage,
not proof that re-evaluated physical densities are unchanged. All outputs are
`reference_frame_candidate`, `qualified=false`, with parent/configuration hashes.

Without `--run`, the tool only prepares candidates; `completed` remains false.
The optional `--carrier-corrector-mode carrier_density` forwards the existing
density-coordinate candidate correction to the direct continuation driver;
the default remains `carrier_only`. Neither mode changes the physical deck or
the full coupled and all-active-row acceptance requirements.
`density_then_qf` uses the first half of the pass budget (rounded up) for density
coordinates, then the remaining passes for QF coordinates, carrying the last
candidate across the transition. It requires at least two passes; eight passes
give at most four density and four QF corrections. An audit or full-validation
failure never grants acceptance; a timeout terminates recovery rather than
switching coordinates. Per-step modes and final validation remain in the ledger.
With `--run --runner ...`, it uses the existing direct driver **serially**:
first fully requalify the optional previous anchor, then fully requalify the
current anchor at the same physical bias, and only then advance requested
targets. A failed previous anchor prevents current/target dispatch. Predictions
use requalified anchors in the same frame and remain candidates themselves.
Every accepted point retains full device, closure and all-active-row auditing.

`--restore-poly-equilibrium` optionally forwards the guarded PolySilicon
candidate restoration after a carrier correction; it requires a positive
`--carrier-corrector-passes` budget. It remains off by default. In addition to
the canonical SG profile, the helper accepts the explicit observed
`ialmob_element_density_grad_qf` profile: that implementation assembles avalanche
only on Si cells. The helper still rejects Poly nodes shared with any other
transport material, floating or conflicting contacts, Poly recombination,
variable Poly mobility and unsupported statistics. All supported cases require
the original full-device and all-row validation after projection. The option
does not remove Poly equations or change tolerances.

`manifest.json` records resource/tool hashes, per-phase results, solver-frame
drain voltage, absolute physical drain voltage and Vds. Referenced physical
resources include the external AverageBox transport CSV; relative paths are
resolved before translation, and changed resources block solver dispatch.
Each phase also retains
the direct driver's per-attempt `reports/adaptive_checkpoint.json`; during an
interruption this may be newer than the completed-phase manifest. Resume by
selecting that independently qualified state and its phase `base.json` as inputs
to a **new** output bundle, with the corresponding physical drain voltage.
Do not overwrite the interrupted bundle or treat a candidate as qualified.

Example (paths are illustrative; no higher-bias targets means anchor verification only):

```powershell
python scripts/run_slot_ldmos_reference_frame.py --base-config SOURCE/base.json --initial-state SOURCE/final_state.h5 --initial-drain-voltage-V 3.125 --target-origin-V 1 --output build/slot_frame_new --runner build-release/vela_example_runner.exe --run
```

Use the UCRT64 runtime PATH and `VELA_LINEAR_SOLVER=umfpack` for the current SLOT
Release case. Default full-Newton budget is 24 iterations; no physical or
acceptance threshold is changed by this wrapper. Returning a qualified state
to another frame can fail numerical qualification and always requires revalidation.

An explicitly selected seed is an initial guess, so its previous bias need not
equal the new target. Mesh, layout, units and potential origin must match; use an
explicit conversion/translation when necessary. The hashed case bundle declares
the source seed and intended model. Loading a seed does not inherit its old
curve qualification.

Electrothermal input/output JSON contains `state_archive.file` and
`state_archive.sha256` instead of large restart arrays. The HDF5 metadata records
which JSON state fields to reconstruct. Nested predictor candidates have their
own references. Production file readers reject inline/archived mixtures and old
inline restart arrays. Scalar JSON reports remain JSON.

Strict electrothermal resume additionally checks the input and mesh snapshots,
sweep and initialization controls, physical-input file hashes, and referenced
accepted/initialization states. A `complete` ledger does not bypass validation.
Resume starts from a saved accepted boundary, not an arbitrary Newton iteration.
Missing lattice temperature is an error, never an implicit 300 K state.

Writes use same-directory temporary files, close HDF5 before replacement, and
check replacement errors. Electrothermal state files are immutable: a writer
refuses to overwrite a file already named by a record. State files are committed
before JSON records and ledgers. An interrupted uncommitted attempt is reported;
it is not silently adopted. Invalid data and replacement failures are tested to
preserve the previous checkpoint. These tests cover process/file failure, not
sudden-power-loss durability.

## Limits and interfaces

Readers validate schema, field names, dimensions, binary type, units and finite
values. Numeric data is bounded to 1 GiB and metadata to 64 KiB. Extra/missing
fields, incomplete split/quantum/density pairs and incorrect modes are rejected.

C++: `StateArchive`, `DDSolutionState`, `ElectrothermalState`, `StateIdentity`.
Python: `state_archive.py`, `electrothermal_state.py`. The migration-only
`migrate_state_seed_to_hdf5.py` imports historical CSV/VDS1 and records parsed
binary64 equality and source/output hashes. It is not a runtime fallback.
Curve/statistics CSV and mesh/doping inputs retain their existing formats.

The SLOT direct-continuation, load-line and Poly-equilibrium drivers use the
same HDF5 contract. Their shared `slot_state_archive.py` checks independent mesh
identity, DD mode, potential origin and nonnegative carrier densities. Secant
and equilibrium outputs remain unqualified candidates with parent hashes;
split QF references and increments are retained separately. Block-step probes
also write HDF5 candidates, never accepted checkpoints. A successful probe
without its expected HDF5 candidate is an interface failure, not a reason to
bisect voltage.
The continuation report records the supplied seed separately; it reports no
qualified bias/state until at least one target passes in that run. A seed's
previous model or accepted-state label does not itself qualify the new campaign.
Timeouts during the coupled solve, row audit or carrier-candidate correction
stop that target without voltage bisection. The report preserves the target and
records the timed-out stage; a missing convergence result is not numerical failure.
After each finished attempt, `reports/adaptive_checkpoint.json` atomically records
the last qualified state and hash, previous predictor state, remaining targets
and last completed attempt. It is flushed once per attempt, not per Newton
iteration. This recovery record does not indicate whether a process is running;
configuration identity and the referenced qualification evidence must be checked
before using it to prepare a new continuation campaign.

SLOT load-line checkpoint schema v2 binds the independently qualified device
file and the circuit confirmation file separately. Numerical confirmation
requires exact decoded field equality; differing container/provenance bytes
alone do not invalidate the state. Both concrete files retain their own hashes.
Legacy CSV checkpoints/v1 load-line ledgers require explicit migration and
requalification, not implicit runtime fallback.

The general codec preserves binary64 values including signed zero. The DD
QF-mode solver writer retains its pre-existing normalization of subnormal values to
zero; population coordinates and their output fields preserve subnormals; migration comparisons distinguish that solver policy from serialization.
HDF5 byte hashes identify concrete files; compare decoded fields when testing
numerical equality between independently written HDF5 files.

See the [migration execution record](validation/templates_ldmos_hdf5_migration_execution_2026-09-23.md)
for tested configurations, repeat timings and platform-specific qualification limits.
