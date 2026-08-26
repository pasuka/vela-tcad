from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "audit_slot_ldmos_delaunay_mesh.py"
SPEC = importlib.util.spec_from_file_location("slot_delaunay", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def mesh_with_opposites(c_y: float, d_y: float) -> dict:
    return {
        "nodes": [
            {"id": 0, "x": 0.0, "y": 0.0},
            {"id": 1, "x": 1.0, "y": 0.0},
            {"id": 2, "x": 0.2, "y": c_y},
            {"id": 3, "x": 0.5, "y": d_y},
        ],
        "regions": [{"id": 0, "name": "Silicon", "material": "Si"}],
        "triangles": [
            {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
            {"id": 1, "node_ids": [1, 0, 3], "region_id": 0},
        ],
    }


class SlotDelaunayAuditTest(unittest.TestCase):
    def run_audit(self, mesh: dict) -> tuple[dict, list[dict]]:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mesh.json"
            path.write_text(json.dumps(mesh), encoding="utf-8")
            summary, edges, _, _ = MODULE.audit(path, (0.5, 0.0), 1.0e-12, 2, 2)
            return summary, edges

    def test_local_obtuse_contribution_can_be_repaired_after_edge_assembly(self) -> None:
        summary, edges = self.run_audit(mesh_with_opposites(0.1, -2.0))
        shared = next(row for row in edges if row["node0"] == 0 and row["node1"] == 1)
        self.assertEqual(shared["local_negative_count"], 1)
        self.assertTrue(shared["delaunay_qualified"])
        self.assertTrue(shared["repaired_by_edge_aggregation"])
        self.assertEqual(summary["counts"]["unresolved_negative_silicon_interior_edges"], 0)

    def test_negative_assembled_dual_face_is_not_delaunay_qualified(self) -> None:
        summary, edges = self.run_audit(mesh_with_opposites(0.1, -0.2))
        shared = next(row for row in edges if row["node0"] == 0 and row["node1"] == 1)
        self.assertLess(shared["signed_dual_length_um"], 0.0)
        self.assertTrue(shared["unresolved_negative_internal"])
        self.assertFalse(summary["decision"]["new_profile_eligible_for_prototyping"])

    def test_boundary_negative_cotangent_is_reported_but_not_called_internal_failure(self) -> None:
        mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": -0.2, "y": 0.1},
            ],
            "regions": [{"id": 0, "name": "Silicon", "material": "Si"}],
            "triangles": [{"id": 0, "node_ids": [0, 1, 2], "region_id": 0}],
        }
        _, edges = self.run_audit(mesh)
        boundary_negative = [
            row for row in edges
            if row["boundary_edge"] and row["local_negative_count"] > 0
        ]
        self.assertTrue(boundary_negative)
        self.assertTrue(all(not row["unresolved_negative_internal"] for row in boundary_negative))


if __name__ == "__main__":
    unittest.main()
