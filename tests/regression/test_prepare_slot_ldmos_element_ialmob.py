"""Physical geometry invariants and isotropic parameter provenance."""
import math
from pathlib import Path
import tempfile
import unittest

from scripts.prepare_slot_ldmos_element_ialmob import native_parameters, triangle_geometry
from scripts.audit_slot_ldmos_element_ialmob import validate_export_mesh


class ElementPreparationTest(unittest.TestCase):
    def test_replay_rejects_equal_sized_but_different_mesh(self):
        mesh={'nodes':[{'id':0,'x':0.,'y':0.},{'id':1,'x':1.,'y':0.},
                       {'id':2,'x':0.,'y':1.}],
              'triangles':[{'id':0,'node_ids':[0,1,2]}]}
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            nodes='id,x_um,y_um\n0,0,0\n1,1,0\n2,0,1\n'
            (root/'nodes.csv').write_text(nodes)
            (root/'elements.csv').write_text('id,node0,node1,node2\n0,0,1,2\n')
            validate_export_mesh(mesh,root)
            (root/'nodes.csv').write_text(nodes.replace('1,1,0','1,1.01,0'))
            with self.assertRaisesRegex(ValueError,'coordinates'): validate_export_mesh(mesh,root)
            (root/'nodes.csv').write_text(nodes)
            (root/'elements.csv').write_text('id,node0,node1,node2\n0,0,2,1\n')
            with self.assertRaisesRegex(ValueError,'topology'): validate_export_mesh(mesh,root)

    def test_geometry_conserves_area_and_scales_without_changing_coefficients(self):
        for points in ([(0,0),(2,0),(.7,1.2)], [(0,0),(2,0),(.1,.2)],
                       [(0,0),(1,0),(0,1)]):
            area, measures, coefficients = triangle_geometry(points)
            self.assertAlmostEqual(sum(measures),area)
            self.assertTrue(all(v >= 0 for v in measures))
            transformed = [(7-3*y,-5+3*x) for x,y in points]
            area2, measures2, coefficients2 = triangle_geometry(transformed)
            self.assertAlmostEqual(area2,9*area)
            for a,b in zip(measures,measures2): self.assertAlmostEqual(9*a,b)
            for a,b in zip(coefficients,coefficients2): self.assertAlmostEqual(a,b)
            # Vertex permutation must permute the corresponding measures.
            _, rotated, _ = triangle_geometry(points[1:]+points[:1])
            for a,b in zip(rotated,measures[1:]+measures[:1]): self.assertAlmostEqual(a,b)

    def test_unnamed_isotropic_parameters_are_shared_without_unit_conversion(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/"original.par"
            p.write_text('IALMob:\n{\n mumin=56,46\n n_ref=9.68e16,2.23e17\n'
                         ' EnormMinimum=0,0\n me_over_m0=1\n mh_over_m0=1.258\n}\n'
                         'IALMob "100":\n{\n mumin=999,999\n}\n')
            params = native_parameters(p)
            e, h = params['electron_parameters_cm'],params['hole_parameters_cm']
            self.assertEqual(e['100'],e['110'])
            self.assertEqual(e['100'],e['111'])
            self.assertEqual(e['100']['mumin'],56)
            self.assertEqual(h['100']['mumin'],46)
            self.assertEqual(e['100']['n_ref'],9.68e16)
            self.assertEqual(h['100']['n_ref'],2.23e17)
            self.assertEqual(e['100']['other_mass'],h['100']['mass'])
            self.assertNotIn('EnormMinimum',e['100'])
            p.write_text('IALMob:\n{\n EnormMinimum=1,0\n}\n')
            with self.assertRaisesRegex(ValueError,'EnormMinimum'): native_parameters(p)

    def test_degenerate_and_nonfinite_geometry_rejected(self):
        for p in ([(0,0),(1,0),(2,0)],[(0,0),(1,0),(0,math.inf)]):
            with self.assertRaises(ValueError): triangle_geometry(p)


if __name__ == '__main__': unittest.main()
