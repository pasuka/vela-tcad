import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bvmethods_qf_paths import qf_drop_paths,integrate_element_vertices


class QfPathsTest(unittest.TestCase):
    def test_common_large_reference_loses_absolute_but_retains_increment_drop(self):
        a=qf_drop_paths([1.,1.],[0.,2.**-60],[[0,1],[1,0]])
        self.assertEqual(a['absolute'][0],0.)
        self.assertEqual(a['anchored'][0],2.**-60)
        self.assertEqual(a['wide'][0],2.**-60)
        self.assertEqual(a['wide'][1],-a['wide'][0])

    def test_anchor_cast_can_lose_small_difference(self):
        a=qf_drop_paths([0.,1.],[1.,2.**-60],[[0,1]])
        self.assertEqual(a['anchored'][0],0.)
        if np.finfo(np.longdouble).nmant>53:
            self.assertEqual(a['wide'][0],2.**-60)

    def test_node_reference_exactly_preserves_representable_drop(self):
        q=np.array([2.**-20,2.**-20+2.**-60])
        a=qf_drop_paths(np.zeros(2),q,[[0,1]])
        b=qf_drop_paths(q,np.zeros(2),[[0,1]])
        for key in a:np.testing.assert_array_equal(a[key],b[key])

    def test_invalid_qf_fields_rejected(self):
        for r,d,e in [([0.],[0.],[[0,1]]),([np.nan],[0.],[[0,0]]),([0.],[0.],[[0.,0.]])]:
            with self.assertRaises(ValueError):qf_drop_paths(r,d,e)

    def test_element_support_conservation_and_order(self):
        args=(np.array([2.,3.,5.,7.]),np.array([.1,.2,.3,.4]),np.array([0,1,0,2]),np.array([4,4,8,8]),np.array([0,1,0,1]))
        a,v=integrate_element_vertices(*args,3)
        self.assertAlmostEqual(a.sum(),sum(args[0]*args[1]))
        np.testing.assert_allclose(v,[.4,.2,.4])
        b,w=integrate_element_vertices(*(x[::-1] for x in args),3)
        np.testing.assert_allclose(a,b);np.testing.assert_allclose(v,w)

    def test_duplicate_or_negative_support_rejected(self):
        with self.assertRaises(ValueError):integrate_element_vertices([1.,2.],[.2,.3],[0,1],[1,1],[0,0],2)
        with self.assertRaises(ValueError):integrate_element_vertices([1.],[-1.],[0],[1],[0],1)

if __name__=='__main__':unittest.main()
