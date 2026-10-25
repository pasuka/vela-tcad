import pathlib
import sys
import unittest
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "scripts"))
from compare_pn2d_staged_zero import field_errors, density_errors


class StagedFieldComparisonTest(unittest.TestCase):
    def test_retains_contact_offset_without_gauge_fitting(self):
        expected = np.array([-.4, 0., .4])
        actual = expected + 7e-6
        metrics = field_errors(actual, expected)
        self.assertAlmostEqual(metrics["max_absolute"], 7e-6, delta=1e-16)
        self.assertAlmostEqual(metrics["rms_absolute"], 7e-6, delta=1e-16)

    def test_density_log_error_covers_both_minority_and_majority(self):
        expected = np.array([1e3, 1e17])
        metrics = density_errors(expected * 10., expected)
        self.assertAlmostEqual(metrics["max_log10_error"], 1.)
        self.assertAlmostEqual(metrics["max_relative"], 9.)

    def test_rejects_invalid_density(self):
        with self.assertRaises(ValueError):
            density_errors(np.array([0.]), np.array([1.]))


if __name__ == "__main__":
    unittest.main()
