import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_g3_rhs_dimension_chain import (
    central_difference,
    drain_current,
    vector_scale,
)


class RhsDimensionChainAuditTest(unittest.TestCase):
    def test_area_factor_invariant_central_difference_has_unit_scale(self) -> None:
        baseline = central_difference({1: 2.0, 2: -1.0}, {1: 6.0, 2: 3.0}, 0.5)
        candidate = central_difference({1: 2.0, 2: -1.0}, {1: 6.0, 2: 3.0}, 0.5)
        result = vector_scale(
            [candidate[node] for node in sorted(candidate)],
            [baseline[node] for node in sorted(baseline)],
        )
        self.assertAlmostEqual(result["least_squares_scale"], 1.0)
        self.assertAlmostEqual(result["cosine"], 1.0)
        self.assertAlmostEqual(result["relative_l2_after_scale"], 0.0)

    def test_drain_current_uses_last_contact_table(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            path = Path(root_text) / "run.log"
            path.write_text(
                " drain 1.000E-01 1.626E-06 -2.490E-14 1.626E-06\n"
                " drain 1.000E-01 3.252E-06 -4.980E-14 3.252E-06\n",
                encoding="utf-8",
            )
            result = drain_current(path)
            self.assertEqual(result["voltage_V"], 0.1)
            self.assertEqual(result["electron_current_A"], 3.252e-6)
            self.assertEqual(result["hole_current_A"], -4.98e-14)


if __name__ == "__main__":
    unittest.main()
