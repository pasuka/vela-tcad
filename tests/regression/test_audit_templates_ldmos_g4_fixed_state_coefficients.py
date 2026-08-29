import unittest

from scripts.audit_templates_ldmos_g4_fixed_state_coefficients import (
    classify_branch,
    variant_config,
)


class TemplatesLdmosG4FixedStateCoefficientAuditTest(unittest.TestCase):
    def test_variants_are_explicit_and_do_not_mutate_baseline(self) -> None:
        baseline = {"solver": {"mobility": {"model": "constant"}}}
        truncated = variant_config(baseline, "truncated_voronoi")

        self.assertNotIn("mesh_geometry", baseline)
        self.assertEqual(
            truncated["mesh_geometry"],
            {
                "node_volume_policy": "barycentric",
                "fallback_negative_cotangent": False,
            },
        )
        self.assertEqual(truncated["solver"]["mobility"]["model"], "constant")

    def test_branch_classifier_prefers_nearest_log_anchor(self) -> None:
        plot_time = classify_branch(
            g4_ratios=[20.0, 11.0, 10.0],
            g3_ratios=[24.5, 11.1, 11.0],
            g4_curve_ratios=[2.2, 1.9, 1.75],
        )
        assembly_real = classify_branch(
            g4_ratios=[2.1, 1.9, 1.8],
            g3_ratios=[24.5, 11.1, 11.0],
            g4_curve_ratios=[2.2, 1.9, 1.75],
        )

        self.assertEqual(plot_time["classification"], "plot_time_supported")
        self.assertEqual(
            assembly_real["classification"], "assembly_real_supported"
        )


if __name__ == "__main__":
    unittest.main()
