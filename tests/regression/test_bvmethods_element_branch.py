"""Properties of a single-cell source observation, independent of branch law."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bvmethods_element_branch import single_element_alpha, observed_vo_transition, observed_vo_alpha_300K


class ElementBranchTests(unittest.TestCase):
    def test_coefficient_is_invariant_to_current_scale(self):
        for scale in [1e-100, 1., 1e100]:
            self.assertAlmostEqual(single_element_alpha([scale]*3, [2.75*scale]*3), 2.75)

    def test_zero_current_does_not_identify_zero_alpha(self):
        self.assertIsNone(single_element_alpha([0.]*3, [0.]*3))
        self.assertEqual(single_element_alpha([3.]*3, [0.]*3), 0.)

    def test_multi_cell_mixture_cannot_be_called_single_element(self):
        with self.assertRaises(ValueError):
            single_element_alpha([1., 2., 3.], [2., 4., 6.])
        with self.assertRaises(ValueError):
            single_element_alpha([1.]*3, [2., 3., 4.])

    def test_invalid_observations_are_rejected(self):
        for unit, source in [([0.]*3, [1.]*3), ([1.]*2, [1.]*2), ([-1.]*3, [1.]*3), ([1.]*3, [float('nan')]*3)]:
            with self.subTest(unit=unit, source=source), self.assertRaises(ValueError):
                single_element_alpha(unit, source)

    def test_transition_is_bounded_monotone_and_has_pure_branches(self):
        import numpy as np
        forces=np.linspace(0, 8e5, 1001)
        w=observed_vo_transition(forces)
        self.assertTrue(np.all(np.diff(w)>=0))
        self.assertTrue(np.all((w>=0)&(w<=1)))
        self.assertEqual(float(observed_vo_transition(3e5)), 0.)
        self.assertEqual(float(observed_vo_transition(5e5)), 1.)
        self.assertAlmostEqual(float(observed_vo_transition(4e5)), .5, delta=2e-10)

    def test_equal_branches_and_zero_force_are_preserved(self):
        import numpy as np
        f=np.array([3.9e5,4e5,4.1e5])
        np.testing.assert_allclose(observed_vo_alpha_300K(f,'e'),7.03e5*np.exp(-1.231e6/f),rtol=3e-16)
        self.assertEqual(float(observed_vo_alpha_300K(0.,'h')), 0.)
        with self.assertRaises(ValueError):observed_vo_transition(-1.)
        with self.assertRaises(ValueError):observed_vo_alpha_300K(4e5,'invalid')


if __name__ == '__main__':
    unittest.main()
