import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_averagebox_node4492 import (
    coefficient_edge,
    local_vela_couple,
    parse_debug_block,
)


class TemplatesLdmosAverageBoxNode4492AuditTest(unittest.TestCase):
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
