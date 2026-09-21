import unittest

from scripts.audit_templates_ldmos_g3_seven_node_interface_ablation import (
    audit,
    audit_endpoint,
)


def carrier(node: int, residual: float, source: float) -> dict[str, str]:
    return {
        "node_id": str(node), "x": str(float(node)), "y": "0",
        "electron_flux": str(residual - source),
        "electron_recombination": str(source), "electron_impact": "0",
        "electron_gauge": "0", "electron_boundary": "0",
        "electron_residual": str(residual),
    }


class TemplatesLdmosG3SevenNodeInterfaceAblationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
                {"id": 3, "x": 1.0, "y": 1.0},
            ],
            "regions": [
                {"id": 0, "name": "silicon", "material": "Silicon"},
                {"id": 1, "name": "oxide", "material": "Oxide"},
            ],
            "triangles": [
                {"id": 10, "region_id": 0, "node_ids": [0, 1, 2]},
                {"id": 11, "region_id": 1, "node_ids": [1, 3, 2]},
            ],
        }
        self.measures = {10: [1.0 / 12.0] * 3, 11: [1.0 / 6.0] * 3}
        self.coefficients = {10: [0.5, 0.0, 0.5], 11: [0.5, 0.0, 0.5]}

    def test_levels_separate_source_volume_and_constrained_pair(self) -> None:
        summary, rows = audit_endpoint(
            bias_v=1.0,
            mesh=self.mesh,
            carrier_rows=[carrier(1, 10.0, 4.0), carrier(2, -6.0, 2.0)],
            measures=self.measures,
            coefficients=self.coefficients,
            transport_materials={"silicon"},
            selected_nodes=(1, 2),
            silicon_potential={1: 0.2, 2: 0.3},
            oxide_potential={1: 0.2, 2: 0.3},
        )
        # Global volume is 1/3 and the Silicon Measure is 1/12: source scale 1/4.
        self.assertAlmostEqual(rows[0]["level_b_region_local_source_residual"], 7.0)
        self.assertAlmostEqual(rows[1]["level_b_region_local_source_residual"], -7.5)
        self.assertEqual(
            rows[0]["level_c_constrained_pair_residual"],
            rows[0]["level_b_region_local_source_residual"],
        )
        self.assertEqual(summary["direct_interface_nodes"], 2)
        self.assertEqual(
            summary["sentaurus_pair_constraint_residual_V"]["maximum_abs"], 0.0
        )

    def test_two_endpoint_report_exposes_growth_and_causal_gates(self) -> None:
        base = {
            "mesh": self.mesh,
            "measures": self.measures,
            "coefficients": self.coefficients,
            "transport_materials": {"silicon"},
            "selected_nodes": (1,),
        }
        report, _ = audit(endpoints=[
            {
                "bias_V": 1.0, "carrier_rows": [carrier(1, 2.0, 0.2)],
                "silicon_potential": {1: 0.1}, "oxide_potential": {1: 0.1},
            },
            {
                "bias_V": 2.0, "carrier_rows": [carrier(1, 6.0, 0.2)],
                "silicon_potential": {1: 0.2}, "oxide_potential": {1: 0.2},
            },
        ], **base)
        self.assertAlmostEqual(report["endpoint_l2_growth_high_over_low"]["level_a"], 3.0)
        self.assertFalse(report["causal_gates"]["level_c_has_independent_frozen_carrier_effect"])
        self.assertFalse(report["contract"]["full_self_consistent_double_node_solve_performed"])

    def test_missing_selected_node_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "lacks selected nodes"):
            audit_endpoint(
                bias_v=1.0, mesh=self.mesh, carrier_rows=[],
                measures=self.measures, coefficients=self.coefficients,
                transport_materials={"silicon"}, selected_nodes=(1,),
                silicon_potential={}, oxide_potential={},
            )


if __name__ == "__main__":
    unittest.main()
