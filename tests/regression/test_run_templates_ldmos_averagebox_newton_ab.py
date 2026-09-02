import json
import unittest
from pathlib import Path

from scripts.run_templates_ldmos_averagebox_newton_ab import (
    FROZEN_HOTSPOTS,
    compare_step_probes,
    with_contact_boundary_reconstruction,
    with_profile,
    without_profile,
)


class TemplatesLdmosAverageBoxNewtonAbTest(unittest.TestCase):
    def test_profile_is_explicit_and_baseline_removes_every_external_field(self) -> None:
        base = {"mesh_geometry": {"node_volume_policy": "barycentric"}}
        candidate = with_profile(base, Path(__file__), 7)
        self.assertEqual(
            candidate["mesh_geometry"]["carrier_transport_couple_profile"],
            "templates_ldmos_external_averagebox",
        )
        self.assertEqual(
            candidate["mesh_geometry"]["external_averagebox_expected_edges"], 7
        )
        baseline = without_profile(candidate)
        self.assertEqual(
            baseline["mesh_geometry"], {"node_volume_policy": "barycentric"}
        )
        self.assertEqual(base, {"mesh_geometry": {"node_volume_policy": "barycentric"}})

    def test_fixed_state_gate_requires_global_and_frozen_hotspot_improvement(self) -> None:
        baseline = []
        candidate = []
        for index, node in enumerate(FROZEN_HOTSPOTS):
            baseline.append({
                "node_id": str(node), "psi_residual": str(index + 1.0),
                "phin_residual": str(10.0 + index),
            })
            candidate.append({
                "node_id": str(node), "psi_residual": str(index + 1.0),
                "phin_residual": str(4.0 + 0.1 * index),
            })
        result = compare_step_probes(baseline, candidate)
        self.assertTrue(result["gate"]["passed"])
        candidate[-1]["phin_residual"] = baseline[-1]["phin_residual"]
        self.assertFalse(compare_step_probes(baseline, candidate)["gate"]["passed"])

    def test_contact_reconstruction_override_is_explicit_and_non_mutating(
        self,
    ) -> None:
        base = {"solver": {"method": "newton"}}
        candidate = with_contact_boundary_reconstruction(
            base, "legacy_node_local")
        self.assertEqual(
            candidate["solver"]["contact_boundary_reconstruction"],
            "legacy_node_local",
        )
        self.assertEqual(base, {"solver": {"method": "newton"}})
        self.assertEqual(with_contact_boundary_reconstruction(base, None), base)

    def test_frozen_diagnostic_contract_records_factorial_and_scope(self) -> None:
        root = Path(__file__).resolve().parents[2]
        path = (
            root / "reference_tcad" / "templates_ldmos_sentaurus2022" /
            "contracts" / "diagnostics" /
            "templates_ldmos_external_averagebox_profile.json"
        )
        contract = json.loads(path.read_text(encoding="utf-8"))
        self.assertFalse(contract["default_enabled"])
        self.assertFalse(
            contract["authorization"]["production_global_default_change"])
        self.assertTrue(contract["authorization"]["phase3_31_point_curve"])
        self.assertEqual(
            contract["authorization"]["curve_gate_status"],
            "fail_maximum_gm_only",
        )
        self.assertEqual(
            len(contract["qualification"]["same_bias_factorial"]), 4)
        qualified = contract["qualification"]["qualified_combination"]
        self.assertEqual(
            qualified["contact_boundary_reconstruction"], "legacy_node_local")
        self.assertLess(
            contract["qualification"]["same_bias_factorial"][-1][
                "log_error_dex"
            ],
            contract["qualification"]["single_bias_stage3_limit_dex"],
        )
        replay = contract["qualification"][
            "same_contract_contact_hfs_frozen_replay"
        ]
        self.assertTrue(replay["passed"])
        self.assertGreater(replay["median_error_improvement_dex"], 1.0)
        feedback = contract["qualification"]["corrected_state_feedback_replay"]
        self.assertLess(feedback["median_feedback_amplification_dex"], 0.01)

    def test_public_schema_names_the_template_private_profile(self) -> None:
        root = Path(__file__).resolve().parents[2]
        schema = json.loads(
            (root / "configs" / "schema" / "vela-simulation.schema.json")
            .read_text(encoding="utf-8")
        )
        geometry = schema["properties"]["mesh_geometry"]["properties"]
        self.assertIn(
            "templates_ldmos_external_averagebox",
            geometry["carrier_transport_couple_profile"]["enum"],
        )
        self.assertEqual(
            geometry["external_averagebox_expected_edges"]["minimum"], 1)
        charge_volume = schema["properties"]["discretization"][
            "properties"
        ]["poisson_charge_volume_policy"]
        self.assertEqual(
            charge_volume["enum"], ["global", "material_local"]
        )


if __name__ == "__main__":
    unittest.main()
