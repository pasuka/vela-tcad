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

The four-equation restart unknowns are potential, the two QF coordinates and
temperature. Density is derived by the physics implementation; the archive does
not invent zero densities to fill a DD-shaped record. Physical QF values and
split reference/increment values are checked for consistency. Keeping the split
arrays separately preserves increments smaller than the ULP of the reference.

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

The general codec preserves binary64 values including signed zero. The DD
solver writer retains its pre-existing normalization of subnormal values to
zero; migration comparisons distinguish that solver policy from serialization.
HDF5 byte hashes identify concrete files; compare decoded fields when testing
numerical equality between independently written HDF5 files.

See the [migration execution record](validation/templates_ldmos_hdf5_migration_execution_2026-09-23.md)
for tested configurations, repeat timings and platform-specific qualification limits.
