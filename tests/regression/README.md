# Regression tests

The engineering example regression suite is driven by `scripts/run_regression.py`
and registered with CTest as the `regression` test. The script copies each
example into a build-local working directory, runs it with `vela_example_runner`,
checks generated outputs, and writes `regression_summary.json`.

The optional Python binding has separate CTest coverage. Configure with
`-DVELA_ENABLE_PYTHON=ON` to build the pybind11 module and register the
`python_api` test, which imports `vela`, loads a mesh, runs Poisson and DC sweep
flows, and verifies generated CSV/VTK outputs.
