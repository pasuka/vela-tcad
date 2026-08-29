import unittest

from scripts.audit_templates_ldmos_phip_srh import p1_gradient


class TemplatesLdmosPhapSrhAuditTest(unittest.TestCase):
    def test_p1_gradient_is_invariant_to_uniform_offset(self) -> None:
        cell = {"node_ids": [0, 1, 2]}
        coordinates = {0: (0.0, 0.0), 1: (2.0, 0.0), 2: (0.0, 3.0)}
        field = {0: 1.0, 1: 5.0, 2: 10.0}
        shifted = {node: value + 0.02 for node, value in field.items()}

        original_gradient = p1_gradient(cell, coordinates, field)
        shifted_gradient = p1_gradient(cell, coordinates, shifted)

        self.assertAlmostEqual(original_gradient[0], shifted_gradient[0])
        self.assertAlmostEqual(original_gradient[1], shifted_gradient[1])


if __name__ == "__main__":
    unittest.main()
