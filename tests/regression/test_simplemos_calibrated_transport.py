"""Checks for conservative allocation and signed error bookkeeping."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import decompose_simplemos_calibrated_transport as d


class AllocationTests(unittest.TestCase):
    def test_shared_interaction_is_allocated_symmetrically(self):
        result=d.allocate({0:0.,1:2.,2:3.,3:9.},('a','b'))
        self.assertAlmostEqual(result['a'],4.)
        self.assertAlmostEqual(result['b'],5.)

    def test_irrelevant_factor_has_zero_contribution(self):
        game={mask:float(2*(mask&1)+3*((mask>>1)&1)) for mask in range(8)}
        result=d.allocate(game,('a','b','unused'))
        self.assertEqual(result['unused'],0.)
        self.assertEqual(sum(result.values()),5.)

    def test_incomplete_game_rejected(self):
        with self.assertRaises(ValueError):d.allocate({0:0.,1:1.},('a','b'))

    def test_cancellation_and_projection_sign(self):
        weights=np.array([[2.,2.,2.]])
        source=np.array([[10.,-9.,0.]])
        self.assertEqual(d.project(weights,source),-2.)
        self.assertEqual(d.project(weights,source)+d.project(weights,-source),0.)

    def test_closed_edge_does_not_change_constant_weight_projection(self):
        self.assertEqual(d.project(np.ones((1,2)),np.array([[3.,-3.]])),0.)


if __name__=='__main__':unittest.main()
