import unittest

from scripts.audit_templates_ldmos_interface_pair_box import (
    audit,
    local_couple,
    region_local_geometry,
)


class TemplatesLdmosInterfacePairBoxAuditTest(unittest.TestCase):
    @staticmethod
    def mesh() -> dict:
        return {
            "nodes": [
                {"id": 0, "x": -1.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 2.0},
                {"id": 3, "x": 0.0, "y": -2.0},
            ],
            "triangles": [
                {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                {"id": 1, "region_id": 1, "node_ids": [1, 0, 3]},
            ],
            "regions": [
                {"id": 0, "name": "silicon", "material": "Si"},
                {"id": 1, "name": "oxide", "material": "SiO2"},
            ],
        }

    def test_region_local_geometry_creates_coincident_potential_slave(self) -> None:
        geometry, pairs = region_local_geometry(self.mesh(), {"si"})
        shared = geometry[(0, 1)]
        self.assertAlmostEqual(
            shared["transport_couple"] / shared["all_couple"], 0.5
        )
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0]["master_dofs"], "psi,phin,phip")
        self.assertEqual(pairs[0]["slave_dofs"], "psi")

    def test_negative_cotangent_uses_production_fallback(self) -> None:
        value, negative = local_couple(
            (-1.0, 0.0), (1.0, 0.0), (0.0, 0.2),
            fallback_negative=True,
        )
        self.assertTrue(negative)
        self.assertAlmostEqual(value, 0.2 / 6.0)

    def test_fixed_state_flux_uses_only_semiconductor_half_box(self) -> None:
        sg_rows = [{
            "edge_id": "0", "node0": "0", "node1": "1",
            "couple_m": "1.5e-6", "electron_flux": "10.0",
        }]
        carrier_rows = []
        for node, sign in ((0, 1.0), (1, -1.0), (2, 0.0), (3, 0.0)):
            carrier_rows.append({
                "node_id": str(node), "x": "0", "y": "0",
                "electron_flux": str(10.0 * sign),
                "electron_flux_abs_sum": "10.0",
                "electron_term_sum": str(10.0 * sign),
                "electron_residual": str(10.0 * sign),
            })
        summary, _, edges, nodes = audit(
            self.mesh(), sg_rows, carrier_rows, {"si"}, {0, 1}
        )
        self.assertAlmostEqual(edges[0]["transport_fraction"], 0.5)
        self.assertAlmostEqual(nodes[0]["candidate_electron_residual"], 5.0)
        self.assertAlmostEqual(
            summary["candidate_over_baseline"]["l2"], 0.5
        )
        self.assertTrue(summary["same_bias_reclose_gate"]["passed"])


if __name__ == "__main__":
    unittest.main()
