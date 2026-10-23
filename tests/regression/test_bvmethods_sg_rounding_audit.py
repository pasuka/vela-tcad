import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from bvmethods_sg_rounding_audit import edge_operand_envelope


class RoundingAuditTests(unittest.TestCase):
    def test_exact_equilibrium_has_no_density_rounding_current(self):
        r=edge_operand_envelope(1e18, 5e20, 0., 0., 1,
                                fixed_q0=True, fixed_q1=True)
        self.assertEqual(r['maximum'], 0.)
        self.assertEqual(r['minimum'], 0.)

    def test_tiny_drive_is_contained_without_macroscopic_change(self):
        r=edge_operand_envelope(5e20, 4.9e20, 0., 1e-14, 1, fixed_q0=True)
        self.assertLessEqual(r['minimum'], r['nominal'])
        self.assertGreaterEqual(r['maximum'], r['nominal'])
        self.assertGreater(r['max_deviation'], 0.)
        self.assertLess(r['max_deviation']/abs(r['nominal']), 1e-12)

    def test_edge_reversal_preserves_current_antisymmetry(self):
        a=edge_operand_envelope(1e17, 2e17, 1e-15, 3e-15, -1)
        b=edge_operand_envelope(2e17, 1e17, 3e-15, 1e-15, -1)
        self.assertAlmostEqual(a['nominal']/-b['nominal'], 1., places=14)
        self.assertAlmostEqual(a['minimum']/-b['maximum'], 1., places=14)

    def test_invalid_operands_fail(self):
        for x in [0., -1., math.nan, math.inf]:
            with self.assertRaises(ValueError):
                edge_operand_envelope(x, 1e18, 0., 1e-14, 1)


if __name__=='__main__':unittest.main()
