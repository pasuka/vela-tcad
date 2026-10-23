import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bvmethods_density_precision import relative_error,summarize


class DensityPrecisionTests(unittest.TestCase):
    def test_small_derivatives_are_not_hidden_by_unit_floor(self):
        self.assertAlmostEqual(relative_error([0.0],[1e-30]),1.0)
        self.assertEqual(relative_error([0.0],[0.0]),0.0)
        self.assertIsNone(relative_error([1e-30],[0.0]))

    def test_norm_does_not_underflow_for_weak_rows(self):
        self.assertAlmostEqual(relative_error([0.0,0.0],[1e-250,1e-250]),1.0)

    def test_loss_reports_diagonal_only_and_preserves_full_check(self):
        rows=[]
        for carrier in [1,2]:
            for block in [0,1,2]:
                rows.append(dict(step='1e-7',carrier=str(carrier),column_block=str(block),row_node='790',column_node='790',
                    J_full='1',J_no_recomb='1',J_difference='0',fd_total='1',fd_term_sum='1',fd_recomb='1e-30'))
        result=summarize(rows)
        self.assertEqual(len(result),6)
        self.assertEqual(result[0]['lost_nonzero_source_diagonals'],[790])
        self.assertEqual(result[0]['full_vs_residual_fd_relative'],0.0)
        self.assertEqual(result[0]['source_difference_vs_direct_fd_relative'],1.0)

    def test_rejects_incomplete_blocks(self):
        with self.assertRaises(ValueError):summarize([{'step':'1e-7','carrier':'1','column_block':'0'}])


if __name__=='__main__':unittest.main()
