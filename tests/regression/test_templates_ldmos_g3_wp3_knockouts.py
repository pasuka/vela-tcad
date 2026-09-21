from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.prepare_templates_ldmos_g3_wp3_knockouts import prepare
from scripts.qualify_templates_ldmos_g3_wp3_knockout import qualify


SOURCE = r'''Electrode {
 { name="drain" Voltage=0.0 }
 { name="gate" Voltage=0.0 }
}
File {
 Grid="n1_fps.tdr"
 Parameters="sdevice.par"
 Output="base.log"
 Current="base.plt"
 Plot="base.tdr"
}
Physics {Fermi}
Physics(Material="Silicon") {
 Mobility(HighFieldSaturation)
 EffectiveIntrinsicDensity(OldSlotboom)
 Recombination(SRH(DopingDependence TempDependence) Auger)
}
Math { Extrapolate Iterations=25 ExitOnFailure }
Solve {
 Coupled { Poisson Electron Hole }
 Quasistationary(Goal { Name="drain" Voltage=0.1 }) { Coupled { Poisson Electron Hole } }
 Quasistationary(Goal { Name="gate" Voltage=1.166666666666667 }) { Coupled { Poisson Electron Hole } }
}
'''


def raw_case(*, ratio: float = 0.05, gate_failure: bool = False, floor_hit: bool = False) -> dict:
    current = 1.0
    baseline = 1.0e-3
    r7_l2_target = baseline * ratio
    flux = 1.0 if not floor_hit else 1.0e-8
    row = r7_l2_target * flux / (7.0 ** 0.5)
    digest = "a" * 64
    return {
        "variant": "k1_hfs_off",
        "bias_V": 1.0,
        "np_log10_errors_dex": {
            "electron": [2.0e-4, 3.0e-4],
            "hole": [4.0e-4, 5.0e-4],
        },
        "contact": {
            "phin_minus_contact_V": [2.0e-10],
            "phip_minus_contact_V": [3.0e-10],
            "neutrality_relative_residual": [2.0e-7],
        },
        "port_current_A_per_um": {
            "sentaurus": current,
            "vela": current * (1.002 if gate_failure else 1.0002),
        },
        "mesh": {
            "vertex_count_equal": True,
            "coordinate_sha256_equal": True,
            "transport_edge_sha256_equal": True,
        },
        "residual": {
            "fixed_seven_rows_A_per_um": [row] * 7,
            "variant_top7_node_ids": [1, 2, 3, 4, 5, 6, 7],
            "variant_top7_rows_A_per_um": [row] * 7,
            "incident_silicon_edge_ids": list(range(40)),
            "incident_silicon_edge_flux_A_per_um": [flux / 40.0] * 40,
        },
        "baseline_normalized_residual": baseline,
        "artifact_sha256": {
            name: digest for name in (
                "deck", "grid", "parameter", "sentaurus_tdr", "imported_state",
                "vela_replay_config", "transport_edges",
            )
        },
    }


class TemplatesLdmosG3Wp3KnockoutTest(unittest.TestCase):
    def test_preparer_emits_three_isolated_two_endpoint_decks_and_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "G3.cmd"
            parameter = root / "sdevice.par"
            grid = root / "n1_fps.tdr"
            source.write_text(SOURCE, encoding="utf-8")
            parameter.write_text("parameter data\n", encoding="utf-8")
            grid.write_bytes(b"sealed synthetic grid")
            manifest = prepare(source, parameter, grid, root / "out")
            self.assertEqual(manifest["status"], "prepared_not_run")
            self.assertEqual(len(manifest["variants"]), 3)
            for record in manifest["variants"]:
                deck = (root / "out" / record["deck"]).read_text(encoding="utf-8")
                self.assertEqual(deck.count("Plot(-Loadable FilePrefix="), 2)
                self.assertIn('Voltage=1.0', deck)
                self.assertIn('Voltage=1.166666666666667', deck)
                self.assertNotIn("IALMob", deck)
                self.assertNotIn("Predictor", deck)
                self.assertEqual(len(record["deck_sha256"]), 64)
                self.assertEqual(record["parameter_copy_sha256"], manifest["parameter_file_sha256"])
                self.assertTrue((root / "out" / record["state_deck_manifest"]).is_file())
                state_manifest = json.loads(
                    (root / "out" / record["state_deck_manifest"]).read_text(encoding="utf-8"))
                self.assertEqual(state_manifest["expected_grid_sha256"], manifest["grid_file_sha256"])
            k1 = (root / "out/k1_hfs_off/IdVg.cmd").read_text(encoding="utf-8")
            k2 = (root / "out/k2_boltzmann/IdVg.cmd").read_text(encoding="utf-8")
            k3 = (root / "out/k3_bgn_off/IdVg.cmd").read_text(encoding="utf-8")
            self.assertNotIn("HighFieldSaturation", k1)
            self.assertNotIn("Mobility()", k1.replace(" ", ""))
            self.assertNotIn("Physics {Fermi}", k2)
            self.assertNotIn("OldSlotboom", k3)
            plan = json.loads((root / "out/execution_plan.json").read_text(encoding="utf-8"))
            self.assertEqual(plan["expected_runs"], 3)
            self.assertEqual(plan["expected_endpoint_states"], 6)

    def test_qualifier_applies_strict_reduction_bands(self) -> None:
        self.assertEqual(qualify(raw_case(ratio=0.05))["verdict"], "discretization_carrier")
        self.assertEqual(qualify(raw_case(ratio=0.75))["verdict"], "family_ruled_out")
        self.assertEqual(qualify(raw_case(ratio=0.2))["verdict"], "uncertain")
        self.assertEqual(qualify(raw_case(ratio=0.1))["verdict"], "uncertain")
        self.assertEqual(qualify(raw_case(ratio=0.5))["verdict"], "uncertain")

    def test_qualification_failure_preempts_reduction_and_floor_is_uncertain(self) -> None:
        failed = qualify(raw_case(ratio=0.05, gate_failure=True))
        self.assertEqual(failed["verdict"], "contract_not_equivalent")
        floored = qualify(raw_case(ratio=0.05, floor_hit=True))
        self.assertTrue(floored["denominator_floor_hit"])
        self.assertEqual(floored["verdict"], "uncertain")

    def test_summary_schema_contains_all_r4_fields_and_hashes(self) -> None:
        schema = json.loads(Path(
            "schemas/vela.templates_ldmos.g3_wp3_knockout_summary.v1.schema.json"
        ).read_text(encoding="utf-8"))
        required = set(schema["required"])
        self.assertTrue({
            "gate_np_p95_dex", "gate_np_max_dex", "gate_contact_qf_max_V",
            "gate_neutrality_max_rel", "gate_port_rel", "gate_mesh_equal",
            "r7_l2", "sum_abs_phi", "denominator_floor_hit",
            "normalized_residual", "ratio_vs_baseline", "hotspot_set", "verdict",
        }.issubset(required))
        self.assertEqual(set(schema["properties"]["artifact_sha256"]["required"]), {
            "deck", "grid", "parameter", "sentaurus_tdr", "imported_state",
            "vela_replay_config", "transport_edges",
        })


if __name__ == "__main__":
    unittest.main()
