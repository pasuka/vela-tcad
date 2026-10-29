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
An independent Linux GCC16 build has positively detected HDF5 state archives,
UMFPACK, SPQR, METIS and the optional MKL PARDISO runtime interface; its full
build and CTest results are still pending. This report does not claim full
Debug CTest completion or a numerical SLOT population-coordinate run.

## Population adaptation status

Not implemented or qualified at this merge boundary. Subsequent changes and
their tests must be recorded separately from the merge baseline. No model
substitution, population floor, protection relaxation or long scan is
authorized by the merge itself. Raw build and simulation evidence remains on
the compute VM; machine-specific records are kept outside version control.
