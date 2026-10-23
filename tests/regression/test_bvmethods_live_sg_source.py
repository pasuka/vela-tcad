import sys, unittest, math
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bvmethods_live_sg_source import (sg_edge_current, representative_current,
                                             van_overstraeten, match_candidates)

class LiveSgSourceTest(unittest.TestCase):
    def test_equilibrium_and_edge_reversal(self):
        self.assertEqual(sg_edge_current(1e10,1e20,0,1e-4,1),0)
        for s in [-1,1]:
            j=sg_edge_current(2e10,3e15,.017,2e-5,s)
            self.assertAlmostEqual(sg_edge_current(3e15,2e10,-.017,2e-5,s)/j,-1,places=13)

    def test_triangle_relabel_and_boundary(self):
        pts=np.array([[0,0],[2e-4,0],[0,1e-4]])
        d=np.array([2e10,3e11,4e10]);q=np.array([.01,.03,-.02]);p=np.array([0,.2,.1])
        j,f=representative_current(pts,d,q,p,0,1)
        perm=[2,0,1]
        jj,ff=representative_current(pts[perm],d[perm],q[perm],p[perm],1,1)
        np.testing.assert_allclose(jj,j,rtol=1e-14);self.assertAlmostEqual(f,ff)
        j,_=representative_current(pts,d,q,p,0,1,boundary_normals=[[1,0]])
        self.assertEqual(j[0],0.)
        contact,_=representative_current(pts,d,q,p,0,1,boundary_normals=[[1,0]],is_contact=True)
        plain,_=representative_current(pts,d,q,p,0,1)
        np.testing.assert_allclose(contact,plain,rtol=0,atol=0)

    def test_duplicate_support_is_not_unique(self):
        cs=[dict(identity=i,current=[1.,2.]) for i in [3,9]]
        m=match_candidates([1.,2.],cs)
        self.assertFalse(m['unique']);self.assertEqual(m['identities'],[3,9])
        self.assertFalse(match_candidates([10.,20.],cs)['unique'])

    def test_live_source_direction_and_field_branch(self):
        self.assertEqual(van_overstraeten(0,'electron'),0.)
        self.assertAlmostEqual(van_overstraeten(4e5,'hole'),6.71e5*math.exp(-1.693e6/4e5))
        self.assertGreater(van_overstraeten(8e5,'electron'),van_overstraeten(4e5,'electron'))
        self.assertLess(van_overstraeten(4e5,'electron',400),van_overstraeten(4e5,'electron',300))
        with self.assertRaises(ValueError):sg_edge_current(-1,2,0,1,1)

if __name__=='__main__':unittest.main()
