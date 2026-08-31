import json
import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_g3_phin_reclose import (
    cosine,
    interpolate_trial_state,
    mesh_node_classes,
    target_step_metrics,
)


class TemplatesLdmosG3PhinRecloseTest(unittest.TestCase):
    def test_interpolate_trial_state_changes_only_carrier_qfs(self) -> None:
        state = [{
            "node_id": 0,
            "psi": 0.2,
            "phin": 0.01,
            "phip": -0.02,
            "electrons_m3": 1.0e20,
            "holes_m3": 2.0e20,
        }]
        step = [{
            "mode": "carrier_only",
            "node_id": "0",
            "trial_phin": "0.05",
            "trial_phip": "-0.04",
        }]
        result = interpolate_trial_state(state, step, 0.25)[0]
        self.assertAlmostEqual(result["psi"], 0.2)
        self.assertAlmostEqual(result["phin"], 0.02)
        self.assertAlmostEqual(result["phip"], -0.025)
        self.assertEqual(result["electrons_m3"], 1.0e20)

    def test_target_step_metrics_detects_aligned_update(self) -> None:
        start = [
            {"node_id": 0, "phin": 0.0},
            {"node_id": 1, "phin": 0.0},
        ]
        target = [
            {"node_id": 0, "phin": 0.1},
            {"node_id": 1, "phin": -0.2},
        ]
        step = [
            {"mode": "carrier_only", "node_id": "0", "delta_phin_V": "0.05"},
            {"mode": "carrier_only", "node_id": "1", "delta_phin_V": "-0.1"},
        ]
        result = target_step_metrics(start, target, step, {0, 1})
        self.assertAlmostEqual(result["phin_update_target_cosine"], 1.0)
        self.assertAlmostEqual(
            result["phin_target_rmse_improvement_fraction"], 0.5
        )

    def test_cosine_handles_zero_direction(self) -> None:
        self.assertEqual(cosine([0.0, 0.0], [1.0, -1.0]), 0.0)

    def test_target_step_metrics_zero_target_has_zero_improvement(self) -> None:
        state = [{"node_id": 0, "phin": 0.0}]
        step = [{
            "mode": "carrier_only",
            "node_id": "0",
            "delta_phin_V": "0.0",
        }]
        result = target_step_metrics(state, state, step, {0})
        self.assertEqual(result["phin_target_rmse_improvement_fraction"], 0.0)

    def test_mesh_node_classes_stratifies_interface_and_contact_ring(self) -> None:
        mesh = {
            "regions": [
                {"id": 0, "material": "Si"},
                {"id": 1, "material": "SiO2"},
            ],
            "triangles": [
                {"region_id": 0, "node_ids": [0, 1, 2]},
                {"region_id": 1, "node_ids": [1, 2, 3]},
            ],
            "contacts": [{"name": "source", "node_ids": [0]}],
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mesh.json"
            path.write_text(json.dumps(mesh), encoding="utf-8")
            classes = mesh_node_classes(path)
        self.assertEqual(classes["interface_nodes"], {1, 2})
        self.assertEqual(classes["interface_one_ring"], {0, 1, 2})
        self.assertEqual(classes["free_nodes"], {1, 2})
        self.assertEqual(classes["contact_one_ring_names"][1], {"source"})


if __name__ == "__main__":
    unittest.main()
