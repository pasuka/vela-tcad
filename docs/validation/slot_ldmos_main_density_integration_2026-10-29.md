# SLOT main integration and population coordinates, 2026-10-29

## Scope and version boundary

Integrate local main `1003a909ce8ffbcf043c30c030c0983c72140cbd` into the SLOT
branch, preserving the pre-integration checkpoint
`8765c3436a8ca1d8ee5ba6818665ddd61962bccf`. Original SLOT physics remains
Boltzmann/IALMob for BVDS and Fermi/PhuMob for IdVd. Integrating main does not
qualify these models for its previously BVMethods-only population runtime.

## Merge verification

All 27 text conflicts were resolved. Semantic integration checks found and
corrected the following issues before committing the merge:

- Accepted load-line trials must reach the existing diagnostic observer before
  state promotion. Observer errors must propagate.
- The SLOT positive-density trial switch and the main boundary-only parameter
  derivative optimization have distinct API arguments.
- Four Fermi perturbation branches must retain scaled dilute-statistics
  evaluation; Joyce--Dixon dispatch remains explicit.
- Named elimination backends must not undergo a second generic elimination
  that removes the retained potential block or changes the reduced cache.
- Restore diagnostic-scale validation and avoid applying a second direction
  limiter to the mapped positive-density trial.
- A first warm target with only one historical state is reported as a constant
  warm start, not as a secant prediction. This changes a reporting assertion,
  not a stopping tolerance.

UCRT64 GCC 16.2 Debug affected targets built successfully. Completed checks:

| Executable or filter | Cases | Assertions | Result |
| --- | ---: | ---: | --- |
| test_newton_solver | 195 | 5581 | pass |
| test_linear_solver | 35 | 675 | pass |
| test_dc_sweep | 124 | 4482 | pass |
| test_carrier_statistics | 33 | 1453 | pass |
| test_bvmethods_population_newton | 11 | 91 | pass |
| test_split_dd_state | 3 | 48 | pass |
| test_split_dd_runtime | 16 | 181 | pass |
| test_phumob | 21 | 4505 | pass |
| test_bvmethods_av_contract | 33 | 990 | pass |
| test_coupled_load_line | 24 | 279 | pass |
| test_element_box_transport [dilute_transport] | 4 | 1347 | pass |

Windows HDF5 test failures under restricted temporary-file replacement were
distinguished from product failures by repeating the same tests with normal
temporary-file permissions. Windows MKL is disabled in this configuration.
An independent Linux GCC16 build positively detected HDF5 state archives,
UMFPACK, SPQR, METIS and the optional MKL PARDISO runtime interface. The first
full merge regression exposed an Eparallel avalanche derivative omission:
the cell electric-field vector and its neighboring-potential derivatives
must also be present without surface-mobility reconstruction. The correction
is covered by the existing residual/Jacobian directional test.

The subsequent population-adaptation build completed all 269 build targets.
Full Linux Release CTest registered 1351 tests: 1346 passed, four native
integration tests were skipped, and `sdevice_package_smoke` failed because
the environment lacks `patchelf`. This is not a clean full-CTest pass, nor a
claim of full Windows Debug CTest completion. The numerical regression tests
passed. Generic compiler-name wrappers and an ignored `build-release/`
temporary directory were needed by two existing regression scripts; these
environment accommodations do not change production solver code.

## Population adaptation status

The merge commit is `e906e7cbabcd9b36c76fe841c4a3f54b926d7c87`. The Eparallel
correction and the following adaptation are later working-tree changes;
they must not be attributed to that merge commit alone.

The opt-in `state_variables: potential_density` path now supports the legacy
Boltzmann/Fermi SLOT operator, in addition to main's BVMethods/Joyce--Dixon
path. This is a coordinate adapter, not a switch to BVMethods physics:

- Internal electron/hole populations are authoritative. Original material,
  mobility, recombination, BGN, avalanche, geometry and statistical inputs
  are retained. Derived QF values are checked against the unchanged relative
  population reconstruction guard.
- Continuity derivatives transform the original physical QF Jacobian into
  potential/population coordinates. The Poisson population columns are
  assembled directly to avoid division of underflowed charge derivatives.
  Statistical inversion uses a relative-accuracy correction on this opt-in
  path; default QF inversion is unchanged.
- Ohmic constraints, terminal derivatives and the external-resistor equation
  use the same coordinate interface. Existing internal-population HDF5
  archives provide authoritative population round trips. Zero-update probes
  evaluate the saved state without writing a corrector update.
- Unsupported split/quantum, thermionic, QF reference/limiter/restoration and
  other coordinate-specific combinations remain explicit errors. The custom
  `slot_population_memory` QF runner is not a population-coordinate runner;
  short validation uses production `dc_sweep`.

Additional UCRT64 checks passed: five legacy-population test cases with 357
assertions, an IALMob geometry/mask/derivative case with 24 assertions, and
the Eparallel directional case with 16 assertions. Tests cover both unit
systems, Boltzmann/Fermi, PhuMob, oxide masks, direct Poisson population
columns, contact/load-line derivatives, inverse accuracy and zero-update
state preservation. Default Newton and DC-sweep regressions also passed.

## Bounded device verification

The resumed two-hour window began on 2026-10-29 at 21:35:58 Beijing time and
ends no later than 23:35:58. Four independent arms compare population and QF
coordinates for IdVd and BVDS. IdVd starts from a matched 6 mV archive and
requests 5/6/7.5 mV; BVDS uses fresh zero-bias Poisson initialization and
requests 0/5/25 mV supply targets with a real 1e12 ohm micrometer resistor.
The device mesh has 13881 nodes; physical port currents are in A/micrometer.
The backend is explicitly MKL PARDISO with generic constraint elimination,
four threads per arm. This does not assert a Sentaurus-matching reduced
matrix dimension.

Both coordinate arms use documented coordinate-compatible numerical settings:
QF-only representation and closure replacement/restoration operations are
disabled, while configured physical closure checks and protection tolerances
remain unchanged. Whole-direction fraction-to-boundary protection is 0.99.
This is not a byte-identical rerun of the historical QF configuration.

The first BVDS submission stopped during configuration parsing: the old
configuration retained `split_contact_bias_low_V` with a drain entry exactly
equal to zero while split state was disabled. A separate retry removes only
this all-zero field, preserves the original failure evidence and shares the
original deadline. Nonzero low parts cannot be discarded by this preparation
step. No production parser guard was relaxed.

The corrected BVDS population arm reached the actual coupled circuit solve,
but rejected its initial zero-supply point during line search. The reported
device equation residual was 2.47094e-12 in the configured equation norm,
while the load-line error was 1.86384e-5 V against the unchanged 1e-6 V
tolerance. No supply point was accepted in that arm. The other arms were
still running at this observation; range completion is not established.
This is an external-circuit convergence failure, not the earlier parser
failure, and does not by itself identify a Jacobian error.

Device-range completion and numerical qualification therefore remain
unproven. Particularly, extreme-population continuity derivatives may lose
information before coordinate transformation; passing the toy/direct-column
tests does not establish full-device convergence in that regime. No model
substitution, population floor, protection relaxation, P2 hard threshold or
long scan was introduced. Raw build and simulation evidence remains on the
compute VM; machine-specific records are kept outside version control.
