import importlib.util
import math
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "audit_templates_ldmos_g3_gm_reclose.py"
SPEC = importlib.util.spec_from_file_location("g3_gm_reclose", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class G3GmRecloseAuditTest(unittest.TestCase):
    def test_vector_comparison_reports_norm_and_direction(self) -> None:
        result = MODULE.vector_comparison([1.0, 2.0], [2.0, 4.0])
        self.assertAlmostEqual(result["actual_over_target_l2"], 0.5)
        self.assertAlmostEqual(result["cosine_similarity"], 1.0)
        self.assertAlmostEqual(result["relative_difference_l2"], 0.5)

    def test_focus_nodes_unions_endpoint_edge_support(self) -> None:
        def document(edges):
            return {"points": [{"edge_feedback": {"top_flux_feedback_edges": edges}}]}

        result = MODULE.focus_nodes([
            document([{"node0": 1, "node1": 2}, {"node0": 2, "node1": 3}]),
            document([{"node0": 3, "node1": 4}]),
        ])
        self.assertEqual(result, [1, 2, 3, 4])

    def test_probe_config_freezes_single_bias_state(self) -> None:
        baseline = {
            "simulation_type": "dc_sweep",
            "sweep": {"bias_points": [0.0, 1.0]},
            "contacts": [{"name": "gate", "bias": 0.0}],
            "solver": {},
        }
        config = MODULE.probe_config(
            baseline, Path("state.csv"), 7.0 / 6.0, Path("probe.csv"),
            "newton_step_probe",
        )
        self.assertNotIn("sweep", config)
        self.assertEqual(config["simulation_type"], "newton_step_probe")
        self.assertTrue(math.isclose(config["contacts"][0]["bias"], 7.0 / 6.0))


if __name__ == "__main__":
    unittest.main()
