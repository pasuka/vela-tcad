from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "audit_slot_ldmos_global_continuity_outputs.py"
SPEC = importlib.util.spec_from_file_location("audit_slot_closure", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SlotClosureAuditTest(unittest.TestCase):
    def write_iv(self, directory: Path, electron_contact: float,
                 electron_source: float, hole_contact: float,
                 hole_source: float, converged: str = "0") -> Path:
        directory.mkdir(parents=True)
        path = directory / "iv.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=[
                "bias_V", "converged", "failure_reason",
                "global_electron_contact_flux", "global_electron_integrated_source",
                "global_hole_contact_flux", "global_hole_integrated_source",
            ])
            writer.writeheader()
            writer.writerow({
                "bias_V": "0.106", "converged": converged,
                "failure_reason": "global_continuity_closure",
                "global_electron_contact_flux": electron_contact,
                "global_electron_integrated_source": electron_source,
                "global_hole_contact_flux": hole_contact,
                "global_hole_integrated_source": hole_source,
            })
        return path

    def test_distinguishes_absolute_candidate_from_true_closure_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = self.write_iv(
                root / "candidate", -1.30e-12, -1.39e-12,
                -1.40e-12, -1.39e-12)
            failed = self.write_iv(
                root / "failed", 9.0e-5, -1.39e-12,
                -1.40e-12, -1.39e-12)
            candidate_result = MODULE.audit_iv(candidate, 1e-11, 1e-12, 0.01)
            failed_result = MODULE.audit_iv(failed, 1e-11, 1e-12, 0.01)
            self.assertEqual(
                candidate_result["classification"],
                "criterion_candidate_rerun_required")
            self.assertTrue(candidate_result["closure_passed"])
            self.assertEqual(failed_result["classification"], "closure_failure")
            self.assertFalse(failed_result["closure_passed"])


if __name__ == "__main__":
    unittest.main()
