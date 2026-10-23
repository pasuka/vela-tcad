import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bvmethods_sg_transfer import (source_rhs, subtraction_resolution,
                                          protected_qf_inverse, block_comparison)


class SgTransferAuditTest(unittest.TestCase):
    def test_source_box_partition_and_units(self):
        # One square micron split among three vertices preserves its source.
        g, nd = 2e20, 1e10
        parts=source_rhs(g,np.array([.2,.3,.5]),nd)
        self.assertAlmostEqual(parts.sum(),200.)
        self.assertAlmostEqual(source_rhs(g,1.,nd),parts.sum())
        np.testing.assert_allclose(source_rhs(g,np.array([0.,1.]),nd),[0.,200.])

    def test_pair_generation_adds_both_triggers(self):
        np.testing.assert_allclose(source_rhs([3.,7.],[2.,5.],11),
                                   source_rhs([1.,2.],[2.,5.],11)+source_rhs([2.,5.],[2.,5.],11))

    def test_invalid_source_inputs_rejected(self):
        for g,m,nd in [(-1.,1.,1.),(1.,-1.,1.),(1.,1.,0.),(np.nan,1.,1.)]:
            with self.assertRaises(ValueError):source_rhs(g,m,nd)

    def test_subtraction_keeps_unresolved_rows(self):
        a=np.array([1.,1.,0.]);b=np.array([np.nextafter(1.,0.),.5,0.])
        d,limit,resolved=subtraction_resolution(a,b)
        self.assertEqual(len(d),3);np.testing.assert_array_equal(resolved,[False,True,False])
        self.assertGreater(d[0],0.);self.assertGreater(limit[0],d[0])

    def test_contact_and_inactive_qf_preserved(self):
        ref=np.array([.1,.2,.3,.4]);req=np.array([8.,9.,np.nan,np.nan])
        out=protected_qf_inverse(ref,req,[1,1,0,0],[1,0,0,1])
        np.testing.assert_array_equal(out,[.1,9.,.3,.4]);np.testing.assert_array_equal(ref,[.1,.2,.3,.4])

    def test_zero_change_identity(self):
        ref=np.array([0.,-.0,5.,np.nextafter(5.,6.)])
        out=protected_qf_inverse(ref,ref,[1,1,1,1],[0,1,0,1])
        self.assertEqual(ref.tobytes(),out.tobytes())

    def test_full_error_does_not_hide_component_cancellation(self):
        z=np.zeros(2);a=np.array([1.,-1.])
        r=block_comparison(z,a,-a,z,z,z,z,z)
        self.assertEqual(r['full_difference_L2'],0.)
        self.assertGreater(r['transport_difference_L2'],1.)
        self.assertEqual(r['internal_closure_relative'],0.)

    def test_bad_shapes_and_nonfinite_rejected(self):
        with self.assertRaises(ValueError):protected_qf_inverse([0.],[0.,1.],[1],[0])
        with self.assertRaises(ValueError):block_comparison(*([[float('nan')]]*8))


if __name__=='__main__':unittest.main()
