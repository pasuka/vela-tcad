import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"scripts"))
from summarize_simplemos_convergence_contrast import paired_change


class PairQualificationTests(unittest.TestCase):
    def test_small_movement_does_not_qualify_failed_states(self):
        result = paired_change(1e-9, 2e-9, True, False)
        self.assertFalse(result["qualified_pair"])

    def test_both_states_must_pass(self):
        self.assertTrue(paired_change(.01, .02, True, True)["qualified_pair"])

    def test_pair_uses_high_minus_low_and_preserves_sign(self):
        self.assertAlmostEqual(paired_change(.03, .01, True, True)["paired_error_change_dex"], -.02)


if __name__ == "__main__":
    unittest.main()
