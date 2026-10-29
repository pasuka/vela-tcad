import unittest
import numpy as np
from scripts.slot_linear_quality import evaluate_coo


class LinearQuality(unittest.TestCase):
    def calc(self, A, b, x, vt=.02585):
        A=np.asarray(A, float);i,j=np.nonzero(A)
        return evaluate_coo(len(b),i,j,A[i,j],b,x,potential_scale_V=vt)

    def test_near_zero_coupled_row_reports_raw_error(self):
        q=self.calc([[1.,0.,0.],[0.,1e-16,-1e-16],[0.,0.,1.]],
                    [1.,0.,0.],[1.,1e-16,0.])
        self.assertTrue(q['accepted']);self.assertEqual(q['raw_bad_count'],1)
        self.assertEqual(q['absolute_allowance_used_count'],1)

    def test_real_weak_row_error_is_rejected(self):
        q=self.calc([[1.,0.,0.],[0.,1e-16,-1e-16],[0.,0.,1.]],
                    [1.,0.,0.],[1.,1e-4,0.])
        self.assertFalse(q['accepted']);self.assertGreater(q['mixed_ratio'],1)

    def test_strong_row_error_is_rejected(self):
        q=self.calc([[2.,1.],[0.,1.]], [3.,1.], [1.1,1.])
        self.assertFalse(q['accepted'])

    def test_row_scaling_invariance(self):
        A=np.array([[2.,1.],[1e-16,-1e-16]])
        b=np.array([3.,0.]);x=np.array([1.,1.+1e-14])
        q=self.calc(A,b,x);scale=np.array([1e-20,1e20])
        r=self.calc(A*scale[:,None],b*scale,x)
        self.assertEqual(q['accepted'],r['accepted'])
        self.assertAlmostEqual(q['mixed_ratio'],r['mixed_ratio'],delta=1e-3)

    def test_coordinate_scale_has_same_physical_allowance(self):
        A=np.array([[2.,1.],[1e-16,-1e-16]]);b=np.array([3.,0.]);x=np.array([1.,1.+1e-14])
        q=self.calc(A,b,x,vt=1.);r=self.calc(A*2,b,x/2,vt=2.)
        self.assertEqual(q['accepted'],r['accepted'])
        self.assertAlmostEqual(q['mixed_ratio'],r['mixed_ratio'],places=10)

    def test_nonfinite_and_inconsistent_zero_row(self):
        self.assertFalse(self.calc([[1.,0.],[0.,1.]], [1.,1.], [np.nan,1.])['accepted'])
        self.assertFalse(self.calc([[1.,0.],[0.,0.]], [1.,1e-30], [1.,0.])['accepted'])

    def test_fixed_constraint_protected_separately(self):
        self.assertTrue(self.calc([[2.,0.],[0.,1.]], [4.,3.], [2.,3.])['accepted'])
        self.assertFalse(self.calc([[2.,0.],[0.,1.]], [4.,3.], [2.+1e-10,3.])['accepted'])

    def test_block_error_cannot_hide_in_absolute_allowance(self):
        A=np.array([[1.,0.,0.,0.],[0.,1.,0.,0.],
                    [0.,0.,1.,-1.],[0.,0.,-1.,1.]])
        i,j=np.nonzero(A)
        q=evaluate_coo(4,i,j,A[i,j],[0.,0.,0.,0.],[0.,0.,1e-20,0.],
                       potential_scale_V=.02585,block_size=2)
        self.assertLess(q['mixed_ratio'],1)
        self.assertEqual(q['block_normwise'][1],1)
        self.assertFalse(q['accepted'])

    def test_duplicate_coordinates_rejected(self):
        with self.assertRaises(ValueError):
            evaluate_coo(1,[0,0],[0,0],[1.,1.],[2.],[1.],potential_scale_V=1.)


if __name__=='__main__':unittest.main()
