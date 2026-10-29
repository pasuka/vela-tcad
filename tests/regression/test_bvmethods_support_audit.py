"""Physical conservation, mapping and source quadrature checks."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from bvmethods_support_audit import edge_divergence,original_vertex_values,element_box_mobility,e2_cell_integrals,e2_nodal_integrals


class SupportAuditTests(unittest.TestCase):
    def test_incidence_conservation_orientation_and_zero_edges(self):
        a=edge_divergence(3,[(0,1),(1,2),(2,0)],[2.,-3.,0.])
        b=edge_divergence(3,[(1,0),(2,1),(0,2)],[-2.,3.,0.])
        np.testing.assert_array_equal(a,b)
        self.assertEqual(sum(a),0.)
        np.testing.assert_array_equal(a,[2.,-5.,3.])

    def test_contact_copies_do_not_contribute_bulk_measure(self):
        values,copies=original_vertex_values([(2,99.),(1,3.),(0,2.)],[[0,0],[1,0]],{0:[0,0],1:[1,0],2:[1,0]})
        np.testing.assert_array_equal(values,[2.,3.]);self.assertEqual(copies,1)

    def test_mapping_rejects_duplicate_missing_and_wrong_coordinates(self):
        for rows,xy in [([(0,1),(0,2)],{0:[0,0]}),([],{0:[0,0]}),([(0,1)],{0:[1,0]})]:
            with self.subTest(rows=rows),self.assertRaises(ValueError):original_vertex_values(rows,[[0,0]],xy)

    def test_box_mean_preserves_constant_mobility_and_transferred_measure(self):
        self.assertEqual(element_box_mobility([0,1,2],[5,5,5],{0:2,1:3,2:9}),5.)
        self.assertAlmostEqual(element_box_mobility([0,1,2],[2,4,8],{0:2,1:3,2:9}),88/14)
        self.assertEqual(element_box_mobility([2,0,1],[8,2,4],{0:2,1:3,2:9}),88/14)

    def test_e2_gauge_permutation_and_geometric_scaling(self):
        points=np.array([[0.,0.],[.1,0.],[0.,.1]]);psi=np.array([0.,2.,3.]);measure=np.array([.002,.001,.002])
        a=e2_cell_integrals(points,psi,measure)
        np.testing.assert_allclose(a,e2_cell_integrals(points,psi+5,measure),rtol=1e-14)
        perm=[2,0,1];np.testing.assert_allclose(a[perm],e2_cell_integrals(points[perm],psi[perm],measure[perm]),rtol=1e-14)
        # Scale geometry and voltage together: field fixed, measure scales as area.
        np.testing.assert_allclose(e2_cell_integrals(points*2,psi*2,measure*4),a*4,rtol=1e-14)
        np.testing.assert_array_equal(e2_cell_integrals(points,[1,1,1],measure),[0,0,0])

    def test_invalid_support_rejected(self):
        with self.assertRaises(ValueError):edge_divergence(2,[(0,2)],[1])
        with self.assertRaises(ValueError):element_box_mobility([0,1,2],[1,-1,1],{0:1,1:1,2:1})
        with self.assertRaises(ValueError):e2_cell_integrals([[0,0],[1,0],[2,0]],[0,1,2],[1,1,1])

    def test_nodal_e2_rotation_measure_and_zero_field(self):
        E=np.array([[3e5,4e5],[0,0.]]);V=np.array([.2,.3]);args=(3.5e21,2.25e7)
        a,J=e2_nodal_integrals(E,V,*args);b,K=e2_nodal_integrals(E[:,::-1],V,*args)
        np.testing.assert_allclose(a,b);np.testing.assert_allclose(J[:,::-1],K)
        np.testing.assert_allclose(e2_nodal_integrals(E,V*4,*args)[0],a*4)
        self.assertEqual(a[1],0);np.testing.assert_array_equal(J[1],[0,0])

    def test_nodal_e2_field_partials(self):
        E=np.array([[3e5,4e5],[2e5,-3e5]]);V=np.array([.2,.3]);direction=np.array([[.1,.2],[-.3,.1]])
        a,J=e2_nodal_integrals(E,V,3.5e21,2.25e7);h=1.
        fd=(e2_nodal_integrals(E+h*direction,V,3.5e21,2.25e7)[0]-e2_nodal_integrals(E-h*direction,V,3.5e21,2.25e7)[0])/(2*h)
        np.testing.assert_allclose(fd,np.sum(J*direction,axis=1),rtol=1e-7)


if __name__=='__main__':unittest.main()
