"""Acceptance boundary checks for the curve extension, without a solver run."""
import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import simplemos_masetti_curve_vela_20260908 as v


class Acceptance(unittest.TestCase):
    def test_both_paths_required_even_when_currents_match(self):
        delta = dict(psi_max_V=0., phin_max_V=0., phip_max_V=0., density_max_relative=0.)
        self.assertTrue(v.dual_qualified(True, True, delta, 0.))
        self.assertFalse(v.dual_qualified(False, True, delta, 0.))
        self.assertFalse(v.dual_qualified(True, False, delta, 0.))

    def test_minority_state_failure_is_not_hidden_by_current_agreement(self):
        delta = dict(psi_max_V=0., phin_max_V=0., phip_max_V=1.01e-6, density_max_relative=0.)
        self.assertFalse(v.dual_qualified(True, True, delta, 0.))
        delta['phip_max_V'] = 0.
        delta['density_max_relative'] = 1.01e-4
        self.assertFalse(v.dual_qualified(True, True, delta, 0.))

    def test_nonfinite_or_large_current_difference_rejects(self):
        delta = dict(psi_max_V=0., phin_max_V=0., phip_max_V=0., density_max_relative=0.)
        for value in (float('nan'), float('inf'), 1.01e-6):
            self.assertFalse(v.dual_qualified(True, True, delta, value))
        for key in delta:
            bad = copy.deepcopy(delta)
            bad[key] = float('nan')
            self.assertFalse(v.dual_qualified(True, True, bad, 0.))

    def test_curve_contains_accepted_anchors_and_both_endpoints_once(self):
        grid = v.n.GRID
        self.assertEqual(len(grid), 51)
        self.assertEqual(len(set(grid)), 51)
        self.assertEqual((grid[0], grid[40], grid[-1]), (0., .8, 1.))


if __name__ == '__main__':
    unittest.main()
