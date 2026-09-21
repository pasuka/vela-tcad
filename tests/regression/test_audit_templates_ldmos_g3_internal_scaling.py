import csv
import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_g3_internal_scaling import (
    errref_formula_check,
    log_summary,
    pearson,
)


def write_field(root: Path, name: str, values: dict[int, float]) -> None:
    path = root / "fields" / f"{name}_region0.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["node_id", "component0"])
        writer.writeheader()
        for node, value in values.items():
            writer.writerow({"node_id": node, "component0": value})


class InternalScalingAuditTest(unittest.TestCase):
    def test_errref_formula_is_recovered_from_exported_error_fields(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            baseline = root / "baseline"
            candidate = root / "candidate"
            density = {1: 0.0, 2: 1.0e10, 3: 1.0e16}
            baseline_error = {1: 2.0, 2: 3.0, 3: 4.0}
            candidate_error = {
                node: baseline_error[node]
                * (abs(value) + 1.0e10) / (abs(value) + 1.0e8)
                for node, value in density.items()
            }
            write_field(baseline, "eDensity", density)
            write_field(baseline, "eDensityError", baseline_error)
            write_field(candidate, "eDensityError", candidate_error)
            result = errref_formula_check(baseline, candidate, 1.0e10, 1.0e8)
            self.assertEqual(result["active_nodes"], 3)
            self.assertLess(result["maximum_relative_error_against_formula"], 1.0e-14)

    def test_log_order_separates_rhs_write_from_first_linear_solve(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            path = Path(root_text) / "run.log"
            path.write_text(
                "RelErrControl (Reference error):\n"
                "  Electron : 1.0000e+10\n"
                "Without diagonal preconditioning\n"
                "Iteration   |Rhs|      factor     |step|     error   #inner\n"
                "    0      7.28e+05\n"
                "writing sample_newton_0_0.tdr (TDR format) ... done.\n"
                "C-norm_equation  max_error vertex coordinate [um] value\n"
                " electron: 3.0e-1 42 (0, 0) 1.0e+12\n"
                "    1      4.72e+09 1.0 1.0 1.0 0 1 1.0\n",
                encoding="utf-8",
            )
            result = log_summary(path)
            self.assertEqual(result["electron_errref_cm3"], 1.0e10)
            self.assertTrue(
                result["iteration0_newtonplot_precedes_cnorm_and_first_linear_solve"]
            )
            self.assertTrue(result["without_diagonal_preconditioning_declared"])

    def test_pearson_reports_no_linear_relation(self) -> None:
        self.assertAlmostEqual(pearson([0.0, 1.0, 2.0], [1.0, 0.0, 1.0]), 0.0)


if __name__ == "__main__":
    unittest.main()
