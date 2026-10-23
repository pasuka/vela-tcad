import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from bvmethods_geometry_field import delaunay_boxes,region_field_stencil
from bvmethods_geometry_field import semiconductor_node_field_stencil
from bvmethods_field_stencil import apply_field_stencil
from bvmethods_field_stencil import e2_potential_derivative
from bvmethods_support_audit import e2_nodal_integrals


class GeometryFieldTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        self.t=[[0,1,2],[0,2,3]];self.r=[3,3]

    def test_affine_field_and_constant_nullspace(self):
        e,_=region_field_stencil(self.p,self.t,self.r,3)
        np.testing.assert_allclose(apply_field_stencil(4,e,self.p@[2.,-3.]+7),[[-2.,3.]]*4,atol=1e-14)
        np.testing.assert_allclose(apply_field_stencil(4,e,[9.]*4),0,atol=1e-14)

    def test_rotation_translation_scaling(self):
        e,_=region_field_stencil(self.p,self.t,self.r,3);q=np.array([[0.,-1.],[1.,0.]])
        rotated,_=region_field_stencil(self.p@q.T*7+[13.,-8.],self.t,self.r,3);psi=[.2,-.3,.5,.7]
        np.testing.assert_allclose(apply_field_stencil(4,rotated,psi),apply_field_stencil(4,e,psi)@q.T/7,atol=1e-14)

    def test_material_support_and_shared_geometry(self):
        p=np.vstack([self.p,[-1.,0.]]);t=self.t+[[0,3,4]]
        e,_=region_field_stencil(p,t,[3,3,0],3)
        self.assertFalse(any(j==4 for i,j,a in e))
        # The other material changes a shared-edge weight, without its potential column.
        e0,_=region_field_stencil(self.p,self.t,self.r,3)
        self.assertGreater(np.linalg.norm(dict(((i,j),a) for i,j,a in e)[0,2]-dict(((i,j),a) for i,j,a in e0)[0,2]),.01)

    def test_obtuse_transfer_preserves_global_coefficient(self):
        p=np.array([[0.,0.],[2.,0.],[1.,.3],[1.,-4.]])
        k,v,n=delaunay_boxes(p,[[0,1,2],[1,0,3]],[3,3])
        self.assertEqual(n,1);self.assertEqual(k[0,0],0);self.assertTrue(np.all(v>=0))
        raw=(.3**2-1)/(4*.3)+(16-1)/16
        self.assertAlmostEqual(k[:,0].sum(),raw)

    def test_non_delaunay_boundary_and_invalid_inputs_reject(self):
        for p,t,r in [([[0,0],[2,0],[1,.1]],[[0,1,2]],[3]),([[0,0],[1,0],[2,0]],[[0,1,2]],[3]),(self.p,[[0,0,1]],[3])]:
            with self.assertRaises(ValueError):delaunay_boxes(p,t,r)
        with self.assertRaises(ValueError):region_field_stencil(self.p,self.t,self.r,8)

    def test_material_boundary_cannot_transfer_negative_part(self):
        p=np.array([[0.,0.],[2.,0.],[1.,.3],[1.,-4.]])
        with self.assertRaises(ValueError):delaunay_boxes(p,[[0,1,2],[1,0,3]],[3,0])
        with self.assertRaises(ValueError):delaunay_boxes(self.p,self.t+[self.t[0]],[3,3,3])

    def test_generated_field_e2_directional_derivative(self):
        e,_=region_field_stencil(self.p*1e-6,self.t,self.r,3)
        psi=np.array([.1,.3,.6,.2]);delta=np.array([.02,-.03,.01,.04]);V=np.ones(4)*.1
        E=apply_field_stencil(4,e,psi);_,p=e2_nodal_integrals(E,V,3.5e21,2.25e7)
        chain=np.zeros(4)
        for i,j,v in e2_potential_derivative(e,p):chain[i]+=v*delta[j]
        h=1e-5
        plus=e2_nodal_integrals(apply_field_stencil(4,e,psi+h*delta),V,3.5e21,2.25e7)[0]
        minus=e2_nodal_integrals(apply_field_stencil(4,e,psi-h*delta),V,3.5e21,2.25e7)[0]
        np.testing.assert_allclose((plus-minus)/(2*h),chain,rtol=1e-7)

    def test_shared_si_vertex_uses_semiconductor_neighbors(self):
        p=np.vstack([self.p,[-1.,0.]]);t=self.t+[[0,3,4]];regions=[3,3,0]
        entries,identity=semiconductor_node_field_stencil(p,t,regions,3)
        self.assertEqual(identity['global_nodes'],5)
        self.assertEqual(identity['semiconductor_nodes'],4)
        self.assertFalse(any(i==0 and j==4 for i,j,a in entries))
        self.assertTrue(any(i==4 and j==0 for i,j,a in entries))
        np.testing.assert_allclose(apply_field_stencil(5,entries,p@[2.,-3.]+7),[[-2.,3.]]*5,atol=1e-14)


if __name__=='__main__':unittest.main()
