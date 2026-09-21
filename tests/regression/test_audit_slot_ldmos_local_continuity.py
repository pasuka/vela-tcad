import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/audit_slot_ldmos_local_continuity.py"
SPEC = importlib.util.spec_from_file_location("local_continuity", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class LocalContinuityAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "terms.csv"

    def row(self, node, residual, active=1):
        row = {"node_id": node, "x": node, "y": 0}
        for carrier in ("electron", "hole"):
            row.update({f"{carrier}_{key}": value for key, value in {
                "continuity_active": active, "residual": residual,
                "flux": residual - 1e-15, "flux_abs_sum": 2e-15,
                "recombination": 1e-15, "impact": 0, "gauge": 0,
                "boundary": 0, "term_sum": residual}.items()})
        return row

    def write(self, rows):
        with self.path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    def test_cancellation_does_not_hide_local_errors_below_source_gate(self):
        self.write([self.row(0, 1e-15), self.row(1, -1e-15), self.row(2, 10, 0)])
        result = MODULE.audit(self.path)
        for carrier in result["carriers"].values():
            self.assertEqual(carrier["active_rows"], 2)
            self.assertEqual(carrier["source_qualified_rows"], 0)
            self.assertEqual(carrier["residual_signed_sum"], 0)
            self.assertEqual(carrier["rows_above_eps"], 2)
            self.assertAlmostEqual(carrier["residual_l1_over_source_l1"], 1)

    def test_nonfinite_or_missing_activity_is_not_silently_excluded(self):
        row = self.row(0, float("nan"))
        self.write([row])
        with self.assertRaises(ValueError):
            MODULE.audit(self.path)
        row = self.row(0, 0)
        del row["electron_continuity_active"]
        self.write([row])
        with self.assertRaises(KeyError):
            MODULE.audit(self.path)

    def test_duplicate_nodes_and_empty_physical_domain_are_invalid(self):
        self.write([self.row(0, 0), self.row(0, 0)])
        with self.assertRaises(ValueError):
            MODULE.audit(self.path)
        self.write([self.row(0, 0, 0)])
        with self.assertRaises(ValueError):
            MODULE.audit(self.path)


if __name__ == "__main__":
    unittest.main()
