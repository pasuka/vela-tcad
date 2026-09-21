import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_averagebox_node4492 import (
    coefficient_edge,
    local_vela_couple,
    parse_debug_block,
    measure_in_tdr_vertex_order,
)
from scripts.audit_templates_ldmos_ialmob_support import project_to_nodes


class TemplatesLdmosAverageBoxNode4492AuditTest(unittest.TestCase):
    def test_measure_places_half_volume_at_right_angle_vertex(self) -> None:
        # Clockwise triangle (0,1), (2,0), (0,0), area=1. The circumcenter
        # decomposition gives vertex volumes 1/4, 1/4, 1/2. Native debug
        # slots in this profile are v0,v2,v1, independently of edge slots.
        mapped = measure_in_tdr_vertex_order([.25,.5,.25])
        self.assertEqual(mapped, [.25,.25,.5])
        self.assertEqual(sum(mapped), 1.)
        with self.assertRaises(ValueError):
            measure_in_tdr_vertex_order([0.,0.,0.])

    def test_plot_projection_uses_inverse_area_and_preserves_constants(self) -> None:
        points={0:(0.,0.),1:(2.,0.),2:(0.,1.),3:(0.,-2.)}
        cells=[{'id':10,'node_ids':[0,1,2]}, {'id':11,'node_ids':[0,3,1]}]
        result=project_to_nodes(points,cells,{10:10.,11:40.})
        self.assertAlmostEqual(result[0],20.)
        self.assertEqual(result[2],10.)
        self.assertEqual(result[3],40.)
        self.assertEqual(set(project_to_nodes(points,cells,{10:7.,11:7.}).values()),{7.})
        with self.assertRaises(ValueError):
            project_to_nodes({**points,2:(1.,0.)},cells,{10:10.,11:40.})

    def test_triangle_coefficient_order_is_tdr_edge_order(self) -> None:
        nodes = [10, 11, 12]
        self.assertEqual(coefficient_edge(nodes, 0), (12, 10, 1))
        self.assertEqual(coefficient_edge(nodes, 1), (11, 12, 0))
        self.assertEqual(coefficient_edge(nodes, 2), (10, 11, 2))

    def test_parse_debug_block_uses_design_element_numbering(self) -> None:
        contents = """
Info { dimension = 2 }
Measure {
  42 7 2 1.0 2.0 3.0
  43 -1 1 9.0 9.0
}
Coefficients {
  42 7 2 4.0 5.0 6.0
}
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "MeasureCoefficients.debug"
            path.write_text(contents, encoding="utf-8")
            self.assertEqual(parse_debug_block(path, "Measure"), {7: [1.0, 2.0, 3.0]})
            self.assertEqual(
                parse_debug_block(path, "Coefficients"), {7: [4.0, 5.0, 6.0]}
            )

    def test_obtuse_triangle_uses_positive_barycentric_fallback(self) -> None:
        couple, policy = local_vela_couple((0.0, 0.0), (1.0, 0.0), (0.5, 0.1))
        self.assertGreater(couple, 0.0)
        self.assertEqual(policy, "positive_barycentric_fallback")


if __name__ == "__main__":
    unittest.main()
