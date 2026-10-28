import sys
import unittest
from pathlib import Path
try:
    import numpy as np
    from scipy import sparse
except ImportError:
    sparse = None
if sparse is not None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/diagnostics"))
    from compare_matrix_support import compare_support


@unittest.skipIf(sparse is None, "optional matrix comparison requires NumPy and SciPy")
class MatrixSupportTests(unittest.TestCase):
    def test_storage_strict_and_common_row_scale(self):
        a = sparse.coo_matrix(([10., 0., 1e-13, 2.], ([0, 0, 1, 1], [0, 1, 0, 1])), shape=(2, 2))
        b = sparse.csr_matrix([[10., 1e-14], [0., 2.]])
        report, positions = compare_support(a, b, (0., 1e-12))
        self.assertEqual((report["left_stored"], report["left_strict_nonzero"]), (4, 3))
        self.assertEqual(report["thresholds"][0]["shared"], 2)
        self.assertEqual(positions[0.]["left_only"], [(1, 0)])
        self.assertEqual(positions[0.]["right_only"], [(0, 1)])
        self.assertEqual(report["thresholds"][1]["left_only"], 0)
        self.assertEqual(report["thresholds"][1]["right_only"], 0)
        self.assertEqual(a.nnz, 4)  # The solver inputs are not pruned.

    def test_row_unit_invariance(self):
        a, b = sparse.csr_matrix([[2., .1], [.3, 4.]]), sparse.csr_matrix([[2., 0.], [.3, 4.]])
        _, expected = compare_support(a, b)
        units = sparse.diags([1e20, 1e-20])
        self.assertEqual(compare_support(units @ a, units @ b)[1], expected)

    def test_rejects_unaligned_or_nonfinite(self):
        with self.assertRaises(ValueError): compare_support(sparse.eye(2), sparse.eye(3))
        with self.assertRaises(ValueError): compare_support(sparse.csr_matrix([[np.nan]]), sparse.eye(1))
        with self.assertRaises(ValueError): compare_support(sparse.eye(2), sparse.eye(2), (-1.,))


if __name__ == "__main__": unittest.main()
