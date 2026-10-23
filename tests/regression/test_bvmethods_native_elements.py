import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from compare_bvmethods_native_elements import population_transform, normalize_system


class NativeElementsTests(unittest.TestCase):
    def test_fixed_population_potential_column_moves_both_qf_columns(self):
        # Different degeneracy slopes require the full chain rule, not log-density alone.
        t = population_transform([0, 0, 0], [0, 1, 2], np.array([[10.], [2.]]),
                                 np.array([[.4], [.8]]), 1, 1.1, 3.)
        np.testing.assert_allclose(t.toarray(), [[1.1, 0, 0], [1.1, -.75, 0], [1.1, 0, 1.875]])
        j = csr_matrix([[2., 3., 5.], [7., 11., 13.], [17., 19., 23.]])
        np.testing.assert_allclose((j @ t).toarray()[:, 0], [11, 34.1, 64.9])

    def test_row_scale_and_rhs_sign_are_both_preserved(self):
        j = csr_matrix([[2., -3.], [5., 7.]])
        rhs = np.array([4., -9.])
        a, b = normalize_system(j, rhs)
        c, d = normalize_system(csr_matrix([[-6., 9.], [25., 35.]]), rhs * [-3., 5.])
        np.testing.assert_allclose(a.toarray(), c.toarray())
        np.testing.assert_allclose(b, d)
        self.assertFalse(np.allclose(b, -d))

    def test_permuted_coordinates_and_invalid_derivative(self):
        t = population_transform([0, 0, 0], [2, 0, 1], np.ones((2, 1)), np.ones((2, 1)), 1, 1, 1)
        np.testing.assert_allclose(t.toarray(), [[0, 1, 0], [0, 1, -1], [1, 1, 0]])
        with self.assertRaises(ValueError):
            population_transform([0], [1], np.zeros((2, 1)), np.ones((2, 1)), 1, 1, 1)
        with self.assertRaises(ValueError):
            normalize_system(csr_matrix([[0.]]), np.array([1.]))


if __name__ == "__main__":
    unittest.main()
