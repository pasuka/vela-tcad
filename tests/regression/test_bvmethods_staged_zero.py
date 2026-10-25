import pathlib
import sys
import tempfile
import unittest
import numpy as np
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / 'scripts'))
from compare_bvmethods_staged_zero import merged_scalar, potential_metrics, carrier_metrics
from compare_bvmethods_staged_operators import poisson_metrics


class StagedBVMethodsTest(unittest.TestCase):
    def test_global_node_ids_and_no_gauge_fit(self):
        m = potential_metrics({4: .002, 12: .302}, {12: .3, 4: 0.})
        self.assertAlmostEqual(m['max_absolute_V'], .002)
        self.assertFalse(m['pass_1mV'])

    def test_minority_population_retained_and_dielectric_excluded(self):
        m = carrier_metrics({4: 1e19, 12: 2., 20: 0.}, {4: 1e19, 12: 1., 20: 0.}, {4, 12})
        self.assertEqual(m['nodes'], 2)
        self.assertEqual(m['worst_node'], 12)
        self.assertFalse(m['pass_1percent'])
        with self.assertRaises(ValueError):
            carrier_metrics({4: 0.}, {4: 1.}, {4})

    def test_conflicting_interface_values_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d)
            for region, value in [(0, .1), (1, .2)]:
                (p/f'Potential_region{region}.csv').write_text(f'node_id,component0\n12,{value}\n')
            with self.assertRaises(ValueError):
                merged_scalar(p, 'Potential')

    def test_poisson_row_and_coordinate_scalings_are_not_mismatch(self):
        a = sparse.csr_matrix([[3., -1.], [-2., 4.]])
        b = np.array([.1, -.2])
        row_scale = sparse.diags([1e-7, 1e5])
        # Candidate x uses twice the reference's volts per packed unit.
        m = poisson_metrics(row_scale@a*2, row_scale@b, a, b, 2., 1.)
        self.assertLess(m['matrix_max_relative_row_L2'], 1e-15)
        self.assertLess(m['rhs_relative_L2'], 1e-15)


if __name__ == '__main__':
    unittest.main()
