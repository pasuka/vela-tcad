import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bvmethods_source_sampling import inverse_current_alpha,reconstruct_edge_mobility,branch_mix_fraction


class SourceSamplingTests(unittest.TestCase):
    def test_inverse_current_cancels_magnitude_and_preserves_zero(self):
        j=np.array([[3.,4.],[-3.,-4.],[0.,0.],[1e-150,0.]])
        a=inverse_current_alpha(j)
        np.testing.assert_allclose(a[:2]*np.hypot(j[:2,0],j[:2,1])/1.602192e-19,1.)
        self.assertEqual(a[2],0.)
        self.assertAlmostEqual(a[3]*1e-150/1.602192e-19,1.)

    def test_constant_mobility_commutes_with_reconstruction(self):
        d=np.array([[1.,0.],[.6,.8]]);j=np.array([2.,-3.])
        np.testing.assert_allclose(reconstruct_edge_mobility(d,j,[7.,7.]),7.*np.linalg.solve(d,j))

    def test_nonuniform_edge_mobility_changes_direction(self):
        j=reconstruct_edge_mobility(np.eye(2),[1.,1.],[2.,8.])
        np.testing.assert_array_equal(j,[2.,8.])
        self.assertNotEqual(j[0],j[1])

    def test_complement_and_fraction(self):
        z=branch_mix_fraction([2.,4.,0.],[2.5,7.,0.],[3.5,5.,0.])
        np.testing.assert_array_equal(z['complement_error'],0.)
        np.testing.assert_array_equal(z['fraction'][:2],[.25,.75])
        self.assertFalse(z['active'][2]);self.assertTrue(np.isnan(z['fraction'][2]))

    def test_branch_response_is_not_clipped(self):
        z=branch_mix_fraction([1.],[3.],[0.])
        self.assertEqual(z['fraction'][0],2.)

    def test_invalid_and_singular_inputs_fail(self):
        with self.assertRaises(ValueError):inverse_current_alpha([[np.nan,0.]])
        with self.assertRaises(ValueError):branch_mix_fraction([1.],[2.,3.],[1.])
        with self.assertRaises(ValueError):reconstruct_edge_mobility(np.eye(2),[1.,2.],[0.,2.])
        with self.assertRaises(np.linalg.LinAlgError):reconstruct_edge_mobility([[1.,0.],[1.,0.]],[1.,2.],[1.,2.])


if __name__=='__main__':unittest.main()
