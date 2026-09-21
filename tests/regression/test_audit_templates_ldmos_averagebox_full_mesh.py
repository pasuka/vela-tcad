import unittest

from scripts.audit_templates_ldmos_averagebox_full_mesh import (
    audit,
    averagebox_geometry,
    ratio,
    transport_profile_rows,
)


class TemplatesLdmosAverageBoxFullMeshAuditTest(unittest.TestCase):
    def test_geometry_aggregates_only_transport_cells(self) -> None:
        mesh = {
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
        couples, measures, volumes, cells, metadata = averagebox_geometry(
            mesh,
            {10: [0.1, 0.2, 0.3], 11: [9.0, 9.0, 9.0]},
            {10: [0.5, 0.0, 0.5], 11: [8.0, 8.0, 8.0]},
            {"silicon"},
        )
        self.assertEqual(metadata["transport_cells"], 1)
        self.assertEqual(metadata["transport_nodes"], 3)
        self.assertEqual(measures, {0: 0.1, 1: 0.2, 2: 0.3})
        self.assertAlmostEqual(couples[(0, 2)], 0.5)
        self.assertAlmostEqual(couples[(0, 1)], 0.5)
        self.assertAlmostEqual(couples[(1, 2)], 0.0)
        self.assertAlmostEqual(volumes[1], 1.0 / 3.0)
        self.assertEqual(len(cells), 3)

    def test_ratio_uses_matching_norm_fields(self) -> None:
        self.assertEqual(
            ratio(
                {"l1": 1.0, "l2": 2.0, "maximum_abs": 3.0},
                {"l1": 2.0, "l2": 4.0, "maximum_abs": 6.0},
            ),
            {"l1": 0.5, "l2": 0.5, "maximum_abs": 0.5},
        )

    def test_runtime_profile_aggregates_cell_local_couples_in_metres(self) -> None:
        rows = transport_profile_rows([
            {"node0": 2, "node1": 1, "averagebox_local_couple_um": 0.25},
            {"node0": 1, "node1": 2, "averagebox_local_couple_um": 0.75},
            {"node0": 0, "node1": 1, "averagebox_local_couple_um": 0.0},
        ])
        self.assertEqual(rows[0], {"node0": 0, "node1": 1, "couple_m": 0.0})
        self.assertEqual(rows[1]["node0"], 1)
        self.assertEqual(rows[1]["node1"], 2)
        self.assertAlmostEqual(rows[1]["couple_m"], 1.0e-6)

    def test_full_replay_gate_uses_all_active_rows_and_frozen_hotspots(self) -> None:
        mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
            ],
            "regions": [
                {"id": 0, "name": "silicon", "material": "Silicon"},
            ],
            "triangles": [
                {"id": 10, "region_id": 0, "node_ids": [0, 1, 2]},
            ],
        }
        sg_rows = [
            {
                "edge_id": "0", "node0": "0", "node1": "1",
                "couple_m": "5e-7", "electron_flux": "10",
                "electron_scaled_flux_per_couple_m": "2e7",
            },
            {
                "edge_id": "1", "node0": "1", "node1": "2",
                "couple_m": "0", "electron_flux": "0",
                "electron_scaled_flux_per_couple_m": "0",
            },
            {
                "edge_id": "2", "node0": "0", "node1": "2",
                "couple_m": "5e-7", "electron_flux": "0",
                "electron_scaled_flux_per_couple_m": "0",
            },
        ]
        carrier_rows = []
        for node, flux, active in ((0, 10.0, True), (1, -10.0, True),
                                   (2, 0.0, False)):
            carrier_rows.append({
                "node_id": str(node), "x": str(mesh["nodes"][node]["x"]),
                "y": str(mesh["nodes"][node]["y"]),
                "electron_flux": str(flux),
                "electron_flux_abs_sum": "10" if active else "0",
                "electron_recombination": "0", "electron_impact": "0",
                "electron_gauge": "0", "electron_boundary": "0",
                "electron_residual": str(flux),
            })

        summary, _, _, nodes = audit(
            mesh, sg_rows, carrier_rows,
            {10: [1.0 / 6.0] * 3},
            {10: [0.25, 0.0, 0.25]},
            {"silicon"}, {0, 1},
        )

        self.assertEqual(summary["row_outcomes"]["active_rows"], 2)
        self.assertAlmostEqual(
            summary["coefficient_only_over_baseline"]["l2"], 0.5
        )
        self.assertFalse(
            summary["hotspot_transfer"]["material_hotspot_transfer_detected"]
        )
        self.assertTrue(summary["external_profile_implementation_gate"]["passed"])
        self.assertAlmostEqual(nodes[0]["coefficient_only_electron_residual"], 5.0)


if __name__ == "__main__":
    unittest.main()
