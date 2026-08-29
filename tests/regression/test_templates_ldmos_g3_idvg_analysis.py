from __future__ import annotations

import unittest

from scripts.analyze_templates_ldmos_g3_idvg import (
    log_intersection,
    max_gm,
    percentile,
)


class TemplatesLdmosG3IdvgAnalysisTest(unittest.TestCase):
    def test_percentile_uses_linear_order_statistic_interpolation(self) -> None:
        self.assertAlmostEqual(percentile([0.0, 1.0, 2.0], 0.95), 1.9)

    def test_fixed_current_intersection_is_linear_in_log_current(self) -> None:
        points = [(0.0, 1.0e-12), (1.0, 1.0e-8)]
        self.assertAlmostEqual(log_intersection(points, 1.0e-10), 0.5)

    def test_max_gm_reports_segment_midpoint(self) -> None:
        result = max_gm([(0.0, 0.0), (1.0, 2.0), (2.0, 3.0)])
        self.assertAlmostEqual(result["value_A_per_um_V"], 2.0)
        self.assertAlmostEqual(result["midpoint_V"], 0.5)


if __name__ == "__main__":
    unittest.main()
