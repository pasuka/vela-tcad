# Regression tests

The engineering example regression suite is driven by `scripts/run_regression.py`
and registered with CTest as the `regression` test. The script copies each
example into a build-local working directory, runs it with `vela_example_runner`,
checks generated outputs, and writes `regression_summary.json`.
