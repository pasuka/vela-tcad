"""Local mobility candidate: dimensional, boundary and derivative properties."""
from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from bvmethods_local_mobility import surface_hfs,triangle_mobility


class LocalMobilityTests(unittest.TestCase):
    def test_limits_monotonicity_and_distance_recovery(self):
        for e in [True,False]:
            self.assertEqual(surface_hfs(400,0,0,1e17,0,e)[0],400)
            value,dn,_=surface_hfs(400,0,0,1e17,0,e)
            fd=(surface_hfs(400,1e-3,0,1e17,0,e)[0]-value)/1e-3
            self.assertAlmostEqual(fd/dn,1,places=4)
            self.assertAlmostEqual(surface_hfs(400,1e5,1e-2,1e17,0,e)[0],400)
            a=surface_hfs(400,1e5,0,1e17,1e5,e)
            self.assertTrue(0<a[0]<400 and a[1]<0 and a[2]<0)

    def test_high_field_velocity_limit(self):
        for e,v in [(True,1.07e7),(False,8.37e6)]:
            self.assertAlmostEqual(surface_hfs(400,0,0,0,1e12,e)[0]*1e12/v,1,places=7)

    def test_scalar_chain_partials_against_centered_difference(self):
        for e in [True,False]:
            a=[400.,1e5,2e-7,1e17,1e5,e];value,dn,df=surface_hfs(*a)
            for index,analytic in [(1,dn),(4,df)]:
                hi=a.copy();lo=a.copy();hi[index]+=.1;lo[index]-=.1
                fd=(surface_hfs(*hi)[0]-surface_hfs(*lo)[0])/.2
                self.assertAlmostEqual(fd/analytic,1,places=7)

    def test_cutoff_is_explicit_and_not_smoothed(self):
        self.assertEqual(surface_hfs(400,0,0,0,np.nextafter(1.,0),True)[2],0)
        self.assertLess(surface_hfs(400,0,0,0,1.,True)[2],0)

    def state(self,contact=False,P=np.eye(2)):
        return [np.array([[0.,0.],[1e-6,0.],[0.,2e-6]]),np.array([0.,.2,.4]),
                np.array([0.,.1,.3]),np.array([0.,1e-6,1e-6]),[400.,300.,200.],[1e17]*3,P,contact,True]

    def test_triangle_live_partials_and_gauge_invariance(self):
        for contact in [True,False]:
            a=self.state(contact,np.diag([0.,1.]));mu,jp,jq=triangle_mobility(*a)
            np.testing.assert_allclose(jp.sum(axis=1),0,atol=1e-10);np.testing.assert_allclose(jq.sum(axis=1),0,atol=1e-10)
            for index,J in [(1,jp),(2,jq)]:
                hi=a.copy();lo=a.copy();direction=np.array([.1,-.2,.3]);h=1e-7
                hi[index]=a[index]+h*direction;lo[index]=a[index]-h*direction
                fd=(triangle_mobility(*hi)[0]-triangle_mobility(*lo)[0])/(2*h)
                np.testing.assert_allclose(fd,J@direction,rtol=2e-7,atol=1e-6)
            shifted=a.copy();shifted[1]=a[1]+2;shifted[2]=a[2]+3
            np.testing.assert_allclose(triangle_mobility(*shifted)[0],mu,rtol=1e-14)

    def test_corner_projection_and_contact_override(self):
        a=self.state(False,np.zeros((2,2)));mu,jp,jq=triangle_mobility(*a)
        np.testing.assert_array_equal(jq,0)
        a[7]=True;contact=triangle_mobility(*a)[0]
        self.assertTrue(np.all(contact<mu))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):surface_hfs(400,1,0,0,-1,True)
        a=self.state();a[6]=np.array([[1,2],[0,1]])
        with self.assertRaises(ValueError):triangle_mobility(*a)


if __name__=='__main__':unittest.main()
