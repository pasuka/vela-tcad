import unittest

from scripts.audit_templates_ldmos_g3_gm_operator import (
    segment_metrics,
    variant_config,
)


class TemplatesLdmosG3GmOperatorTest(unittest.TestCase):
    def test_variants_are_explicit_and_non_mutating(self) -> None:
        base = {
            "solver": {
                "mobility": {
                    "model": "constant_field",
                    "high_field_driving_force": "quasi_fermi_gradient",
                    "contact_electric_field_fallback": True,
                    "contact_electric_field_fallback_scope": "contact_node_cell",
                    "contact_electric_field_fallback_mode":
                        "cell_gradient_magnitude",
                }
            },
            "mesh_geometry": {
                "carrier_transport_couple_profile":
                    "templates_ldmos_external_averagebox",
                "external_averagebox_couples_file": "couples.csv",
                "external_averagebox_expected_edges": 10,
            },
        }

        mesh = variant_config(base, "mesh_default")
        contact_off = variant_config(base, "contact_hfs_off")
        global_off = variant_config(base, "global_hfs_off")

        self.assertNotIn(
            "carrier_transport_couple_profile", mesh["mesh_geometry"]
        )
        self.assertEqual(mesh["solver"]["mobility"]["model"], "constant_field")
        self.assertFalse(
            contact_off["solver"]["mobility"][
                "contact_electric_field_fallback"
            ]
        )
        self.assertEqual(global_off["solver"]["mobility"], {"model": "constant"})
        self.assertTrue(
            base["solver"]["mobility"]["contact_electric_field_fallback"]
        )

    def test_segment_metric_uses_exact_two_point_slope(self) -> None:
        result = segment_metrics({1.0: 2.0, 7.0 / 6.0: 3.0}, 7.5)
        self.assertAlmostEqual(result["gm_A_per_um_V"], 6.0)
        self.assertAlmostEqual(result["ratio_to_sentaurus_gm"], 0.8)
        self.assertAlmostEqual(result["relative_error_to_sentaurus_gm"], 0.2)
        self.assertAlmostEqual(result["endpoint_current_growth_ratio"], 1.5)


if __name__ == "__main__":
    unittest.main()
