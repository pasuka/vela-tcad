import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bvmethods_qf_transport_split import divergence, qf_components, residual_partition, VT


class QfTransportSplitTests(unittest.TestCase):
    def test_contact_flux_is_retained_and_internal_flux_cancels(self):
        flux=divergence(3,[0,1],[1,2],[2.,3.])
        np.testing.assert_array_equal(flux,[2.,1.,-3.])
        self.assertEqual(flux.sum(),0.)
        # Removing the contact endpoint first would destroy the terminal budget.
        self.assertEqual(flux[1:].sum(),-2.)
        np.testing.assert_array_equal(divergence(3,[1,2],[0,1],[-2.,-3.]),flux)

    def test_unequal_DOS_and_voltage_conventions_close_exact_QF_split(self):
        psi=np.array([.3,.7]);ni=np.array([1e10,2e10]);dos=np.array([3e19,4e19])
        eta=np.array([-2.,3.]); ne=np.array([-1.9,3.1]);band=np.array([-.4,-.6]);q=np.array([.02,-.01])
        for sign in [-1,1]:
            parts=qf_components(psi,ni,dos,band,q,eta,ne,VT*.99999,sign)
            expected=psi+sign*VT*(eta-np.log(ni/dos))-q
            np.testing.assert_allclose(parts.sum(axis=0),expected,rtol=1e-14,atol=1e-16)
            self.assertTrue(np.all(parts[1]!=0))
            self.assertTrue(np.all(parts[2]!=0))

    def test_vector_partition_keeps_cancellation_and_all_tail_rows(self):
        a=np.array([100.,1e-20]);b=np.array([2.,3e-20]);c=np.array([3.,-1e-20]);d=np.array([4.,2e-20]);n=np.array([1.,0.])
        parts=residual_partition(a,b,c,d,n)
        np.testing.assert_allclose(parts.sum(axis=0),a-n,rtol=1e-14,atol=0.)
        self.assertEqual(parts.shape,(4,2))
        self.assertGreater(sum(np.linalg.norm(x) for x in parts),np.linalg.norm(a-n))

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError): divergence(2,[0],[2],[1.])
        with self.assertRaises(ValueError): divergence(2,[0],[1],[float('nan')])
        with self.assertRaises(ValueError): residual_partition([1.],[1.,2.],[1.],[1.],[1.])
        with self.assertRaises(ValueError): qf_components(0.,0.,1.,0.,0.,0.,0.,VT,-1)


if __name__=='__main__': unittest.main()
