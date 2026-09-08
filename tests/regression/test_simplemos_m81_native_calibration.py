"""Guard the sign and amplitude conventions used to qualify native IFM."""
import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from analyze_simplemos_m81_fixed_charge_fd import compare
from run_simplemos_m81_output_amendment import rows


class NativeCalibrationTests(unittest.TestCase):
    def test_quadratic_dc_error_cancels_in_central_difference(self):
        # I(q) = 2 + 3q + 0.2q^2, q = density / 1e13.
        r = compare(3.55, .55, 2., 3., 5e12, 2.)
        self.assertAlmostEqual(r["relative_error"], 0.)
        self.assertTrue(r["positive_perturbation_sign_pass"])
        self.assertTrue(r["negative_perturbation_sign_pass"])
        self.assertGreater(r["even_nonlinear_fraction"], 0.)

    def test_circuit_sign_must_not_be_fitted_away(self):
        r = compare(2.1, 1.9, 2., -.1, 1e13, 2.)
        self.assertFalse(r["same_sign"])
        self.assertAlmostEqual(r["relative_error"], 2.)

    def test_baseline_drift_remains_visible(self):
        r = compare(2.1, 1.9, 2., .1, 1e13, 1.99)
        self.assertAlmostEqual(r["signal_to_zero_control_drift"], 10.)

    def test_zero_signal_rejected(self):
        with self.assertRaises(ValueError):
            compare(2., 2., 2., .1, 1e13, 2.)

    def test_duplicate_native_datasets_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "native.plt"
            path.write_text('Info { datasets = ["frequency" "frequency"] } Data { 0 0 }')
            with self.assertRaises(ValueError):
                rows(path)


if __name__ == "__main__":
    unittest.main()
