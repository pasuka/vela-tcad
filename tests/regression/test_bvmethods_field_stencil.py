import unittest
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from bvmethods_field_stencil import apply_field_stencil, e2_potential_derivative
from bvmethods_support_audit import e2_nodal_integrals


class FieldStencilTests(unittest.TestCase):
    def setUp(self):
        # Exact gradient of a right triangle, with a constant offset nullspace.
        self.entries = [(i,j,a) for i in range(3) for j,a in enumerate(([1,1],[-1,0],[0,-1]))]

    def test_affine_reproduction_and_gauge(self):
        np.testing.assert_allclose(apply_field_stencil(3,self.entries,[3,5,7]),[[-2,-4]]*3)
        np.testing.assert_array_equal(apply_field_stencil(3,self.entries,[7,7,7]),np.zeros((3,2)))

    def test_linearity_and_rotation(self):
        a=np.array([.4,.1,-.3]); b=np.array([-.2,.3,.7])
        np.testing.assert_allclose(apply_field_stencil(3,self.entries,2*a-b),2*apply_field_stencil(3,self.entries,a)-apply_field_stencil(3,self.entries,b))
        Q=np.array([[0,-1],[1,0]])
        rotated=[(i,j,Q@v) for i,j,v in self.entries]
        np.testing.assert_allclose(apply_field_stencil(3,rotated,a),apply_field_stencil(3,self.entries,a)@Q.T)

    def test_e2_potential_chain_rule(self):
        E=np.array([[4e5,1e5],[2e5,-3e5],[-5e5,1e5]])
        V=np.array([.1,.2,.3]); delta=np.array([.13,-.27,.41])
        entries=[(i,j,np.array(a)*1e4) for i,j,a in self.entries]
        _,partials=e2_nodal_integrals(E,V,3.5e21,2.25e7)
        expected=np.zeros(3)
        for i,j,x in e2_potential_derivative(entries,partials):expected[i]+=x*delta[j]
        dE=apply_field_stencil(3,entries,delta);h=1e-4
        plus=e2_nodal_integrals(E+h*dE,V,3.5e21,2.25e7)[0]
        minus=e2_nodal_integrals(E-h*dE,V,3.5e21,2.25e7)[0]
        np.testing.assert_allclose((plus-minus)/(2*h),expected,rtol=2e-8)

    def test_zero_field_and_volume(self):
        _,p=e2_nodal_integrals(np.zeros((3,2)),[1,0,1],3.5e21,2.25e7)
        self.assertTrue(all(x==0 for i,j,x in e2_potential_derivative(self.entries,p)))

    def test_invalid_entries_fail(self):
        for entries in [[(0,3,[1,2])],[(0,0,[1,np.nan])],[(0,0,[1,2])]*2]:
            with self.assertRaises(ValueError):apply_field_stencil(3,entries,[0,0,0])
        with self.assertRaises(ValueError):e2_potential_derivative([],[[1,2,3]])


if __name__=='__main__':unittest.main()
