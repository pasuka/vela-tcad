from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "audit_templates_ldmos_g3_wp3_t5_hfs_targeted.py"
SPEC = importlib.util.spec_from_file_location("t5_hfs", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class T5HfsTargetedTests(unittest.TestCase):
    def test_magnitude_average_is_area_weighted_and_not_vector_cancelled(self) -> None:
        cells = [
            {
                "area2": 1.0,
                "electron_gradient_norm_V_per_um": 2.0,
                "hole_gradient_norm_V_per_um": 3.0,
                "touches_contact": False,
            },
            {
                "area2": 3.0,
                "electron_gradient_norm_V_per_um": 4.0,
                "hole_gradient_norm_V_per_um": 5.0,
                "touches_contact": True,
            },
        ]
        field, contact = MODULE.magnitude_average_drive(cells, "electron")
        self.assertAlmostEqual(field, 3.5e6)
        self.assertTrue(contact)

    def test_field_limit_matches_frozen_hfs_formula(self) -> None:
        mu0 = 0.1417
        vsat = 107000.0
        beta = 1.109
        self.assertEqual(MODULE.field_limited_mobility(mu0, 0.0, vsat, beta), mu0)
        expected = mu0 / (1.0 + (mu0 * 1.0e6 / vsat) ** beta) ** (1.0 / beta)
        self.assertAlmostEqual(
            MODULE.field_limited_mobility(mu0, 1.0e6, vsat, beta), expected
        )

    def test_transport_reconstruction_applies_carrier_dirichlet_rows(self) -> None:
        rows = [
            {"node0": "0", "node1": "1", "electron_flux": "2.0"},
            {"node0": "1", "node1": "2", "electron_flux": "3.0"},
        ]
        sums, absolute = MODULE.reconstruct_transport(rows, 3, {2})
        self.assertEqual(sums, {0: 2.0, 1: 1.0, 2: 0.0})
        self.assertEqual(absolute, {0: 2.0, 1: 5.0, 2: 0.0})

    def test_point_gate_cannot_be_hidden_by_aggregate_improvement(self) -> None:
        baseline = {node: 1.0 for node in MODULE.FIXED_NODES}
        candidate = {node: 0.1 for node in MODULE.FIXED_NODES}
        candidate[MODULE.FIXED_NODES[0]] = 0.6
        for node in range(100, 114):
            baseline[node] = 10.0
            candidate[node] = 9.0
        sets = {
            "silicon_nodes": set(baseline),
            "interface_band": set(baseline),
        }
        result = MODULE.evaluate_gate(baseline, candidate, sets)
        self.assertLess(result["fixed_seven"]["ratio"]["l2"], 0.5)
        self.assertFalse(result["fixed_seven"]["passes_half_gate"])
        self.assertFalse(result["passes_t5_gate"])


if __name__ == "__main__":
    unittest.main()
