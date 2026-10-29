import sys
from pathlib import Path
import unittest
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from constraint_elimination import eliminate_constraints


class ConstraintEliminationTests(unittest.TestCase):
    def test_nonzero_prescribed_contact_increment_and_circuit_border(self):
        # x1 is constrained to a nonzero affine function of circuit voltage x2.
        a = sparse.csc_matrix([[4., 2., -1.], [0., 2., -3.], [1., 5., 7.]])
        b = np.array([3., 4., -2.])
        reduced, rhs, kept, recover = eliminate_constraints(a, b, [1])
        x = recover(spsolve(reduced, rhs))
        np.testing.assert_allclose(a @ x, b, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(x, spsolve(a, b), rtol=1e-14, atol=1e-14)
        np.testing.assert_array_equal(kept, [0, 2])
        # Simply dropping the constrained row/column fails this physical load.
        dropped = np.zeros(3)
        dropped[kept] = spsolve(a[kept, :][:, kept], b[kept])
        self.assertGreater(np.linalg.norm(a @ dropped - b), 1.)

    def test_rejects_coupled_or_singular_constraint_block(self):
        for a in [sparse.csc_matrix([[1., 1.], [1., 1.]]), sparse.csc_matrix((2, 2))]:
            with self.assertRaises(ValueError):
                eliminate_constraints(a, np.zeros(2), [0, 1])

    def test_rejects_duplicate_constraint_ids(self):
        with self.assertRaises(ValueError):
            eliminate_constraints(sparse.eye(2), np.ones(2), [0, 0])


if __name__ == '__main__':
    unittest.main()
