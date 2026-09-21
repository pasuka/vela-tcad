import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "audit_templates_ldmos_g3_hfs_gradient_ab.py"
SPEC = importlib.util.spec_from_file_location("hfs_gradient_ab", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AB = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AB)


def rows(values):
    return [
        {"node_id": str(node), "electron_residual": str(value)}
        for node, value in values.items()
    ]


class HfsGradientAbTest(unittest.TestCase):
    def setUp(self):
        self.nodes = set(AB.FIXED_NODES) | {8000, 8001}
        self.node_sets = {
            "silicon": self.nodes,
            "free_silicon": self.nodes,
            "interface_band": set(AB.FIXED_NODES) | {8000},
        }

    def test_pass_requires_half_gate_and_no_migration(self):
        baseline = {node: float(index + 1) for index, node in enumerate(sorted(self.nodes))}
        candidate = {node: 0.4 * value for node, value in baseline.items()}
        result = AB.analyze_endpoint(rows(baseline), rows(candidate), self.node_sets)
        self.assertTrue(result["fixed_seven"]["passes_half_gate"])
        self.assertTrue(result["passes_migration_guard"])
        self.assertTrue(result["passes_t3_gate"])

    def test_hotspot_migration_fails_guard(self):
        baseline = {node: 10.0 for node in self.nodes}
        candidate = {node: 4.0 for node in self.nodes}
        candidate[8001] = 20.0
        result = AB.analyze_endpoint(rows(baseline), rows(candidate), self.node_sets)
        self.assertTrue(result["fixed_seven"]["passes_half_gate"])
        self.assertFalse(result["passes_migration_guard"])
        self.assertFalse(result["passes_t3_gate"])


if __name__ == "__main__":
    unittest.main()
