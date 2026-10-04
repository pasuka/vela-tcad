"""Check calibration preserves thresholds and the physical intervention."""
from pathlib import Path
import sys
import unittest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m79b_numerical_calibration as m


class CalibrationEvidence(unittest.TestCase):
    def test_no_retrospective_threshold_relaxation(self):
        if not m.CONTRACT.exists():
            self.skipTest("run after contract freeze")
        self.assertEqual(m.old.read_json(m.CONTRACT)["thresholds"],
                         m.old.read_json(m.old.CONTRACT)["thresholds"])

    def test_original_failure_record_is_retained(self):
        report = m.old.verify()
        self.assertFalse(report["m80_released"])
        self.assertIn("pass_duality", report["failed_gates"])

    def test_qualification_preserves_source_and_current(self):
        statuses = list(m.LOCAL.glob("n*/vd_*/*/adjoint.status.json"))
        if not statuses:
            self.skipTest("run after calibration")
        for path in statuses:
            oldpath = m.old.LOCAL / path.relative_to(m.LOCAL)
            if not oldpath.exists():
                continue
            calibrated, original = m.old.read_json(path), m.old.read_json(oldpath)
            self.assertEqual(calibrated["current_A_per_um"], original["current_A_per_um"])
            self.assertEqual(calibrated["electron_volume_response"]["source_norm"],
                             original["electron_volume_response"]["source_norm"])
            self.assertEqual(calibrated["electron_volume_response"]["continuity_source_norm"], 0)
            self.assertEqual(calibrated["electron_volume_response"]["parameter_direct_current_A_per_um"], 0)
            self.assertEqual(calibrated["electron_volume_response"]["new_nonlinear_solves"], 0)


if __name__ == "__main__":
    unittest.main()
