from __future__ import annotations

import unittest

from scripts.prepare_slot_ldmos_jvp_audit import build_document


class SlotLdmosJvpAuditPreparationTest(unittest.TestCase):
    def test_builds_node_local_step_sweep_and_sets_inner_bias(self) -> None:
        base = {
            "simulation_type": "dc_sweep",
            "output_csv": "iv.csv",
            "contacts": [
                {"name": "source", "bias": 0.0},
                {"name": "drain", "bias": 0.0},
            ],
            "solver": {
                "impact_ionization": {"source_jacobian": "frozen"}
            },
            "sweep": {"contact": "drain"},
        }
        document = build_document(
            base,
            state_file="state.csv",
            output_csv="jvp.csv",
            hotspot_node=10236,
            drain_bias_V=15.720950570922257,
        )
        self.assertEqual(document["simulation_type"], "newton_jvp_probe")
        self.assertNotIn("sweep", document)
        self.assertEqual(
            document["solver"]["impact_ionization"]["source_jacobian"],
            "local_ad",
        )
        drain = next(
            contact for contact in document["contacts"]
            if contact["name"] == "drain"
        )
        self.assertEqual(drain["bias"], 15.720950570922257)
        self.assertEqual(len(document["directions"]), 9)
        self.assertEqual(document["directions"][0]["node_ids"], [10236])
        self.assertEqual(
            {row["amplitude_V"] for row in document["directions"]},
            {1.0e-4, 1.0e-6, 1.0e-8},
        )

    def test_supports_avalanche_free_low_voltage_deck(self) -> None:
        base = {
            "simulation_type": "dc_sweep",
            "contacts": [
                {"name": "source", "bias": 0.0},
                {"name": "drain", "bias": 0.0},
            ],
            "solver": {"mobility": {"model": "masetti_field_lombardi"}},
            "sweep": {"contact": "drain"},
        }
        document = build_document(
            base,
            state_file="state.csv",
            output_csv="jvp.csv",
            hotspot_node=6445,
            drain_bias_V=1.0,
            disable_impact=True,
        )
        self.assertNotIn("impact_ionization", document["solver"])


if __name__ == "__main__":
    unittest.main()
