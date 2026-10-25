"""Properties of the retained contact-row operation used for native comparison."""
import pathlib
import sys
import unittest

import numpy as np
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "scripts"))
from compare_pn2d_zero_operator import canonical_contacts, equation_rhs_difference


class ContactCanonicalizationTest(unittest.TestCase):
    def test_preserves_nonzero_constraint_rhs_and_circuit_direction(self):
        # Two prescribed potential rows, one free device row, one circuit row.
        matrix = np.array([[2., 0., 0., -2.], [3., 8., -4., 7.],
                           [0., 0., 4., 0.], [5., -3., 6., 11.]])
        rhs = np.array([.12, -.4, -.08, .3])
        canonical, transformed_rhs = canonical_contacts(
            sparse.csc_matrix(matrix), rhs, np.array([0, 2]))
        np.testing.assert_allclose(np.linalg.solve(canonical.toarray(), transformed_rhs),
                                   np.linalg.solve(matrix, rhs), rtol=2e-14, atol=2e-15)
        np.testing.assert_allclose(canonical.toarray()[[1, 3]][:, [0, 2]], 0.)

    def test_rejects_coupled_contact_block(self):
        matrix = sparse.csc_matrix([[2., 1., 3.], [0., 4., 5.], [1., 2., 6.]])
        with self.assertRaisesRegex(ValueError, "diagonal"):
            canonical_contacts(matrix, np.ones(3), np.array([0, 1]))

    def test_rhs_unit_conversion_is_independent_of_candidate_jacobian(self):
        physical_rhs = np.array([2., -7., 3.])
        vela_units = np.array([4., .2, 8.])
        native_units = np.array([.5, 3., 2.])
        # Different approximate Jacobians must not create a false RHS error.
        native_diagonal = native_units * np.array([9., 2., 4.])
        vela_diagonal = vela_units * np.array([5., 7., 11.])
        vb, nb = physical_rhs * vela_units, physical_rhs * native_units
        self.assertGreater(np.linalg.norm(vb / vela_diagonal - nb / native_diagonal), .1)
        np.testing.assert_allclose(equation_rhs_difference(
            vb, nb, native_diagonal, native_units / vela_units), 0., atol=1e-15)


if __name__ == "__main__":
    unittest.main()
