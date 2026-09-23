from __future__ import annotations

import copy
import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.prepare_slot_ldmos_vela_ialmob_ablation import (
    build_probe_case,
    count_interface_edges,
    normalized_physics,
    prepare,
)


BASE = {
    "solver": {
        "mobility": {
            "model": "masetti_field",
            "doping_concentration_basis": "total_impurity",
            "high_field_driving_force": "quasi_fermi_gradient",
        },
        "impact_ionization": {"coupling_mode": "self_consistent"},
        "handoff": {"gummel_max_iter": 50},
    },
    "output_csv": "outputs/stages/05_avalanche_on_60v/iv.csv",
    "sweep": {
        "bias_points": [10.0, 60.0],
        "start": 10.0,
        "stop": 60.0,
        "initial_state_file": "outputs/stages/04/final_state.h5",
        "external_circuit": {"initial_inner_voltage_V": 0.01},
        "boundary_control": {
            "resume": True,
            "evaluation_csv": (
                "outputs/stages/05_avalanche_on_60v/evaluations.csv"
            ),
        },
        "write_state_file": "outputs/stages/05_avalanche_on_60v/final_state.h5",
    },
}


class SlotLdmosVelaIALMobAblationTest(unittest.TestCase):
    def test_counts_exact_region_pair_edges(self) -> None:
        mesh = {
            "regions": [
                {"id": 0, "name": "Silicon_1"},
                {"id": 1, "name": "Oxide_1"},
            ],
            "triangles": [
                {"region_id": 0, "node_ids": [0, 1, 2]},
                {"region_id": 1, "node_ids": [1, 0, 3]},
            ],
        }
        self.assertEqual(count_interface_edges(mesh, "Silicon_1", "Oxide_1"), 1)

    def test_pair_changes_only_mobility_and_output_namespace(self) -> None:
        off = build_probe_case(copy.deepcopy(BASE), "ialmob_off", 0.733)
        on = build_probe_case(copy.deepcopy(BASE), "ialmob_on", 0.733)
        self.assertEqual(off["solver"]["mobility"]["model"], "masetti_field")
        self.assertEqual(
            on["solver"]["mobility"]["model"], "masetti_field_lombardi"
        )
        self.assertEqual(
            on["solver"]["mobility"]["surface"]["surface_interface"],
            ["Silicon_1", "Oxide_1"],
        )
        self.assertEqual(normalized_physics(off), normalized_physics(on))
        self.assertEqual(off["solver"]["handoff"]["gummel_max_iter"], 0)
        self.assertEqual(on["solver"]["handoff"]["gummel_max_iter"], 0)

    def test_resume_is_explicit_and_shared(self) -> None:
        off = build_probe_case(copy.deepcopy(BASE), "ialmob_off", 0.733, True)
        on = build_probe_case(copy.deepcopy(BASE), "ialmob_on", 0.733, True)
        self.assertTrue(off["sweep"]["boundary_control"]["resume"])
        self.assertTrue(on["sweep"]["boundary_control"]["resume"])

    def test_prepare_creates_runtime_output_directories(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bundle = Path(directory)
            mesh = {
                "regions": [
                    {"id": 0, "name": "Silicon_1"},
                    {"id": 1, "name": "Oxide_1"},
                ],
                "triangles": [
                    {"region_id": 0, "node_ids": [0, 1, 2]},
                    {"region_id": 1, "node_ids": [1, 0, 3]},
                ],
            }
            (bundle / "mesh.json").write_text(json.dumps(mesh), encoding="utf-8")
            (bundle / "simulation_05_avalanche_on_60v.json").write_text(
                json.dumps(BASE), encoding="utf-8"
            )
            iv = bundle / "outputs/stages/05_avalanche_on_60v/iv.csv"
            iv.parent.mkdir(parents=True)
            with iv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["converged", "outer_voltage_V", "inner_voltage_V"],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "converged": "1",
                        "outer_voltage_V": "60",
                        "inner_voltage_V": "0.733",
                    }
                )
            prepare(bundle)
            for case in ("ialmob_off", "ialmob_on"):
                root = bundle / "outputs/ialmob_ablation/probe_60v" / case
                self.assertTrue((root / "boundary_control_checkpoints").is_dir())
                self.assertTrue((root / "states").is_dir())
                self.assertTrue((root / "vtk").is_dir())


if __name__ == "__main__":
    unittest.main()
