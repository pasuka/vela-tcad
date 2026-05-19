#!/usr/bin/env python3
"""Unit tests for regression-runner policy helpers."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import stat
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("run_regression", REPO / "scripts" / "run_regression.py")
assert SPEC is not None and SPEC.loader is not None
run_regression = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_regression)


FIELDNAMES = [
    "mode",
    "bias_contact",
    "bias_V",
    "current_contact",
    "current_electron",
    "current_hole",
    "current_total",
    "converged",
    "iterations",
    "step_diagnostics",
]


def base_row(**updates: str) -> dict[str, str]:
    row = {
        "mode": "iv",
        "bias_contact": "anode",
        "bias_V": "0",
        "current_contact": "cathode",
        "current_electron": "0",
        "current_hole": "0",
        "current_total": "0",
        "converged": "1",
        "iterations": "1",
        "step_diagnostics": "attempted_step=1;accepted_step=1;retry_count=0",
    }
    row.update(updates)
    return row


def write_example(tmp: Path, cfg: dict[str, object], rows: list[dict[str, str]], fieldnames: list[str] | None = None) -> Path:
    example_dir = tmp / "example"
    (example_dir / "outputs").mkdir(parents=True)
    (example_dir / "simulation.json").write_text(json.dumps(cfg) + "\n")
    with (example_dir / "outputs" / "sweep.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames or FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return example_dir


class RegressionRunnerPolicies(unittest.TestCase):
    def test_ldmos_fieldplate_trend_accepts_ratio_within_limit(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "bv_reverse", "start": 0, "stop": 0.1, "step": 0.1},
            "regression": {
                "ldmos_fieldplate_trend": {
                    "baseline_config": "simulation_bv.json",
                    "baseline_csv": "outputs/baseline.csv",
                    "max_field_ratio_limit": 1.2,
                }
            },
        }
        fieldnames = FIELDNAMES + ["max_electric_field_V_per_m"]
        rows = [
            base_row(mode="bv_reverse", max_electric_field_V_per_m="8"),
            base_row(mode="bv_reverse", bias_V="0.1", current_total="2e-12", max_electric_field_V_per_m="10"),
        ]
        baseline_cfg = {
            "output_csv": "outputs/baseline.csv",
            "sweep": {"mode": "bv_reverse", "start": 0, "stop": 0.1, "step": 0.1},
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows, fieldnames)
            (example_dir / "simulation_bv.json").write_text(json.dumps(baseline_cfg) + "\n")
            with (example_dir / "fake_runner.sh").open("w") as f:
                f.write("#!/bin/sh\ncp outputs/sweep.csv outputs/baseline.csv\n")
            (example_dir / "fake_runner.sh").chmod(stat.S_IRWXU)
            result = run_regression.check_ldmos_fieldplate_trend(example_dir, example_dir / "fake_runner.sh")
            self.assertLessEqual(result["max_field_ratio"], 1.2)

    def test_ldmos_fieldplate_trend_rejects_excessive_ratio(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "bv_reverse", "start": 0, "stop": 0.1, "step": 0.1},
            "regression": {
                "ldmos_fieldplate_trend": {
                    "baseline_config": "simulation_bv.json",
                    "baseline_csv": "outputs/baseline.csv",
                    "max_field_ratio_limit": 1.1,
                }
            },
        }
        fieldnames = FIELDNAMES + ["max_electric_field_V_per_m"]
        rows = [base_row(mode="bv_reverse", max_electric_field_V_per_m="10")]
        baseline_cfg = {"output_csv": "outputs/baseline.csv", "sweep": {"mode": "bv_reverse", "start": 0, "stop": 0, "step": 1}}
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows, fieldnames)
            (example_dir / "simulation_bv.json").write_text(json.dumps(baseline_cfg) + "\n")
            with (example_dir / "fake_runner.sh").open("w") as f:
                f.write("#!/bin/sh\ncat > outputs/baseline.csv <<'CSV'\nmode,bias_contact,bias_V,current_contact,current_electron,current_hole,current_total,converged,iterations,step_diagnostics,max_electric_field_V_per_m\nbv_reverse,anode,0,cathode,0,0,0,1,1,attempted_step=1;accepted_step=1;retry_count=0,5\nCSV\n")
            (example_dir / "fake_runner.sh").chmod(stat.S_IRWXU)
            with self.assertRaisesRegex(AssertionError, "max field ratio"):
                run_regression.check_ldmos_fieldplate_trend(example_dir, example_dir / "fake_runner.sh")
    def test_monotone_current_tolerance_scales_to_observed_currents(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "iv", "start": 0, "stop": 1, "step": 1},
            "regression": {"dc_sweep": {"expected_rows": 2, "require_monotone_abs_current": True}},
        }
        rows = [base_row(current_total="8e-13"), base_row(bias_V="1", current_total="0")]
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows)
            with self.assertRaisesRegex(AssertionError, "abs\\(current_total\\).+decreased"):
                run_regression.check_dc_sweep_regression(example_dir)

    def test_bv_max_field_monotone_failure_names_field_and_row(self) -> None:
        bv_columns = FIELDNAMES + [
            "max_electric_field_V_per_m",
            "current_jump_ratio",
            "breakdown_detected",
            "breakdown_voltage",
            "criterion",
            "last_stable_bias",
            "failed_bias",
            "failure_reason",
        ]
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "bv_reverse", "start": 0, "stop": 1, "step": 1},
            "regression": {
                "dc_sweep": {
                    "expected_rows": 2,
                    "require_monotone_max_field": True,
                    "max_field_monotone_abs_tolerance": 0.0,
                    "max_field_monotone_rel_tolerance": 0.0,
                }
            },
        }
        rows = [
            base_row(
                mode="bv_reverse",
                max_electric_field_V_per_m="10",
                current_jump_ratio="0",
                breakdown_detected="0",
                breakdown_voltage="0",
                criterion="",
                last_stable_bias="0",
                failed_bias="0",
                failure_reason="",
            ),
            base_row(
                mode="bv_reverse",
                bias_V="1",
                max_electric_field_V_per_m="9",
                current_jump_ratio="1",
                breakdown_detected="0",
                breakdown_voltage="0",
                criterion="",
                last_stable_bias="0",
                failed_bias="0",
                failure_reason="",
            ),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows, bv_columns)
            with self.assertRaisesRegex(AssertionError, "row 2 field=max_electric_field_V_per_m"):
                run_regression.check_dc_sweep_regression(example_dir)

    def test_declared_not_converged_allows_nonconverged_rows(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "iv", "start": 0, "stop": 0, "step": 1},
            "regression": {"declared_converged": False, "dc_sweep": {"expected_rows": 1}},
        }
        rows = [base_row(converged="0")]
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows)
            run_regression.check_csv_converged(example_dir)
            self.assertEqual(run_regression.check_dc_sweep_regression(example_dir)["converged_rows"], 0)

    def test_final_bv_nonconvergence_exception_is_shared_by_csv_and_dc_checks(self) -> None:
        bv_columns = FIELDNAMES + [
            "max_electric_field_V_per_m",
            "current_jump_ratio",
            "breakdown_detected",
            "breakdown_voltage",
            "criterion",
            "last_stable_bias",
            "failed_bias",
            "failure_reason",
        ]
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "bv_reverse", "start": 0, "stop": -1, "step": -1},
            "regression": {"dc_sweep": {"expected_rows": 2, "allow_nonconverged_final_bv_point": True}},
        }
        rows = [
            base_row(
                mode="bv_reverse",
                max_electric_field_V_per_m="1",
                current_jump_ratio="0",
                breakdown_detected="0",
                breakdown_voltage="",
                criterion="",
                last_stable_bias="",
                failed_bias="",
                failure_reason="",
            ),
            base_row(
                mode="bv_reverse",
                bias_V="-1",
                converged="0",
                max_electric_field_V_per_m="2",
                current_jump_ratio="1",
                breakdown_detected="1",
                breakdown_voltage="-1",
                criterion="last_stable_before_nonconvergence",
                last_stable_bias="0",
                failed_bias="-1",
                failure_reason="solver failed",
            ),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows, bv_columns)
            run_regression.check_csv_converged(example_dir)
            self.assertEqual(run_regression.check_dc_sweep_regression(example_dir)["converged_rows"], 1)

    def test_ldmos_iv_rejects_empty_csv(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "iv", "start": 0, "stop": 0.1, "step": 0.1},
            "regression": {"ldmos_iv": {}},
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, [])
            with self.assertRaisesRegex(AssertionError, "LDMOS DD-IV CSV contains no rows"):
                run_regression.check_ldmos_iv_trend(example_dir)

    def test_ldmos_iv_rejects_non_monotone_abs_current(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "iv", "start": 0, "stop": 0.2, "step": 0.1},
            "regression": {
                "ldmos_iv": {
                    "current_monotone_abs_tolerance": 1e-30,
                    "current_monotone_rel_tolerance": 0.0,
                }
            },
        }
        rows = [
            base_row(current_total="0"),
            base_row(bias_V="0.1", current_total="2e-12"),
            base_row(bias_V="0.2", current_total="1e-12"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows)
            with self.assertRaisesRegex(AssertionError, "LDMOS \\|Id\\|-Vd trend is not monotone"):
                run_regression.check_ldmos_iv_trend(example_dir)

    def test_ldmos_iv_uses_relative_monotonicity_tolerance(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "iv", "start": 0, "stop": 0.2, "step": 0.1},
            "regression": {
                "ldmos_iv": {
                    "drain_current_sign": 1.0,
                    "current_monotone_abs_tolerance": 1e-30,
                    "current_monotone_rel_tolerance": 0.1,
                }
            },
        }
        rows = [
            base_row(current_total="0"),
            base_row(bias_V="0.1", current_total="1.0e-9"),
            base_row(bias_V="0.2", current_total="9.5e-10"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows)
            result = run_regression.check_ldmos_iv_trend(example_dir)
            self.assertEqual(result["abs_currents"], [0.0, 1.0e-9, 9.5e-10])

    def test_ldmos_iv_rejects_wrong_drain_current_sign(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "iv", "start": 0, "stop": 0.1, "step": 0.1},
            "regression": {"ldmos_iv": {"drain_current_sign": 1.0}},
        }
        rows = [base_row(current_total="0"), base_row(bias_V="0.1", current_total="-1e-12")]
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows)
            with self.assertRaisesRegex(AssertionError, "does not match expected polarity"):
                run_regression.check_ldmos_iv_trend(example_dir)

    def test_ldmos_iv_rejects_missing_required_columns(self) -> None:
        cfg = {
            "output_csv": "outputs/sweep.csv",
            "sweep": {"mode": "iv", "start": 0, "stop": 0.1, "step": 0.1},
            "regression": {"ldmos_iv": {}},
        }
        fieldnames = [name for name in FIELDNAMES if name != "current_total"]
        rows = [base_row(), base_row(bias_V="0.1")]
        for row in rows:
            row.pop("current_total")
        with tempfile.TemporaryDirectory() as tmpdir:
            example_dir = write_example(Path(tmpdir), cfg, rows, fieldnames)
            with self.assertRaisesRegex(AssertionError, "missing column 'current_total'"):
                run_regression.check_ldmos_iv_trend(example_dir)

    def test_run_example_tolerates_configured_nonzero_runner_exit(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            repo = tmp / "repo"
            source = repo / "examples" / "dev"
            source.mkdir(parents=True)
            cfg = {
                "output_csv": "outputs/sweep.csv",
                "sweep": {"mode": "iv", "start": 0, "stop": 0, "step": 1},
                "regression": {"declared_converged": False, "dc_sweep": {"expected_rows": 1}},
            }
            (source / "simulation.json").write_text(json.dumps(cfg) + "\n")
            runner = tmp / "runner.py"
            runner.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "Path('outputs').mkdir(exist_ok=True)\n"
                "Path('outputs/sweep.csv').write_text("
                "'mode,bias_contact,bias_V,current_contact,current_electron,current_hole,current_total,converged,iterations,step_diagnostics\\n'"
                "+ 'iv,anode,0,cathode,0,0,0,0,1,attempted_step=1;accepted_step=1;retry_count=0\\n')\n"
                "raise SystemExit(1)\n"
            )
            runner.chmod(runner.stat().st_mode | stat.S_IXUSR)
            result = run_regression.run_example(
                runner,
                repo,
                tmp / "work",
                {
                    "name": "dev",
                    "config": Path("examples/dev/simulation.json"),
                    "expected": [Path("outputs/sweep.csv")],
                    "checks": ["csv_converged", "dc_sweep_regression"],
                },
            )
            self.assertTrue(result["passed"], result)
            self.assertTrue(result["runner_nonzero_exit_allowed"])


if __name__ == "__main__":
    unittest.main()
