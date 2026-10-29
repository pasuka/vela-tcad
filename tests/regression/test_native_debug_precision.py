import io
import pathlib
import sys
import unittest

import numpy as np
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "scripts"))
from native_debug_precision import rounding_radius, variable_radii, propagate_state_rounding


class DebugPrecisionTest(unittest.TestCase):
    def test_preserves_printed_trailing_zeros_and_exponent(self):
        self.assertEqual(rounding_radius("-1.56139203128890E+01"), 5e-14)
        self.assertEqual(rounding_radius("0.00000000000000E+00"), 5e-15)
        self.assertEqual(rounding_radius("1.200E-03"), 5e-7)
        self.assertNotEqual(rounding_radius("1.2"), rounding_radius("1.200"))

    def test_iteration_and_variable_block_selection(self):
        stream = io.StringIO("Newton loop 0\nVAR:Variable potential:\nVAR:#0 = 1.000E+01\n"
                             "Newton loop 0\nVAR:Variable potential:\nVAR:#0 = 2.00E-01\n"
                             "VAR:Variable density:\nVAR:#7 = 1.000E+04\n")
        self.assertEqual(variable_radii(stream, 1), {(0, 0): .0005, (1, 7): 5.})

    def test_linear_propagation_encloses_all_box_corners(self):
        matrix = sparse.csr_matrix([[2., -7., 3.], [-4., 0., -8.]])
        radius = np.array([.02, .1, .003])
        bound = propagate_state_rounding(matrix, radius)
        for bits in range(8):
            perturbation = radius * [1 if bits & (1 << j) else -1 for j in range(3)]
            self.assertTrue(np.all(abs(matrix @ perturbation) <= bound + 2e-16))

    def test_rejects_invalid_radii(self):
        with self.assertRaises(ValueError):
            propagate_state_rounding(sparse.eye(2), [1., -1.])


if __name__ == "__main__":
    unittest.main()
