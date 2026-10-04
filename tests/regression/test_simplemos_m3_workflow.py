from __future__ import annotations

import importlib.util
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "simplemos_m3_workflow",
    ROOT / "scripts" / "run_simplemos_m3_workflow.py",
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
CONTRACT = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_sdevice_validation_contract_v1.json"
)


class SimpleMosM3WorkflowTest(unittest.TestCase):
    def base(self, root: Path) -> Path:
        for name in ("mesh.json", "materials.json"):
            (root / name).write_text("{}\n", encoding="utf-8")
        (root / "doping.csv").write_text(
            "node_id,donors,acceptors\n", encoding="utf-8")
        path = root / "base.json"
        path.write_text(json.dumps({
            "simulation_type": "dc_sweep",
            "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv",
            "materials_file": "materials.json",
            "contacts": [
                {"name": "source", "bias": 0.0},
                {"name": "drain", "bias": 0.0},
                {"name": "gate", "bias": 0.0},
                {"name": "substrate", "bias": 0.0},
            ],
            "solver": {"method": "newton", "max_iter": 20},
            "sweep": {"mode": "iv"},
        }) + "\n", encoding="utf-8")
        return path

    def test_normalized_controls_convert_to_physical_voltage(self) -> None:
        drain = MODULE.physical_step_control({
            "InitialStep": 0.1,
            "Increment": 1.5,
            "MinStep": 1.0e-5,
            "MaxStep": 1.0,
        }, 0.0, 0.05)
        self.assertAlmostEqual(drain["initial_step_V"], 0.005)
        self.assertAlmostEqual(drain["min_step_V"], 5.0e-7)
        self.assertAlmostEqual(drain["max_step_V"], 0.05)

        gate = MODULE.physical_step_control({
            "InitialStep": 0.01,
            "Increment": 1.5,
            "MinStep": 1.0e-5,
            "MaxStep": 0.05,
        }, 0.0, 2.5)
        self.assertAlmostEqual(gate["initial_step_V"], 0.025)
        self.assertAlmostEqual(gate["min_step_V"], 2.5e-5)
        self.assertAlmostEqual(gate["max_step_V"], 0.125)

    def test_invalid_normalized_control_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "normalized interval"):
            MODULE.physical_step_control({
                "InitialStep": 1.1,
                "Increment": 1.5,
                "MinStep": 1.0e-5,
                "MaxStep": 1.0,
            }, 0.0, 1.0)

    def test_materialized_branches_are_strict_three_stage_chains(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m3_") as temporary:
            root = Path(temporary)
            manifest = MODULE.materialize(
                self.base(root), CONTRACT, root / "run", [0.05, 1.0])

            self.assertEqual("materialized", manifest["status"])
            self.assertEqual(6, len(manifest["stages"]))
            self.assertTrue(manifest["rules"]["previous_accepted_state_only"])
            for offset in (0, 3):
                equilibrium, drain, gate = manifest["stages"][offset:offset + 3]
                self.assertIsNone(equilibrium["predecessor"])
                self.assertIsNone(equilibrium["initial_state_file"])
                self.assertEqual(equilibrium["id"], drain["predecessor"])
                self.assertEqual(
                    equilibrium["final_state_file"], drain["initial_state_file"])
                self.assertEqual(drain["id"], gate["predecessor"])
                self.assertEqual(
                    drain["final_state_file"], gate["initial_state_file"])

                equilibrium_cfg = MODULE.read_json(Path(equilibrium["config"]))
                drain_cfg = MODULE.read_json(Path(drain["config"]))
                gate_cfg = MODULE.read_json(Path(gate["config"]))
                self.assertEqual(
                    {"mode": "poisson_block"},
                    equilibrium_cfg["sweep"]["initialization"])
                self.assertNotIn("initial_state_file", equilibrium_cfg["sweep"])
                self.assertEqual(
                    equilibrium["final_state_file"],
                    drain_cfg["sweep"]["initial_state_file"])
                self.assertEqual(
                    drain["final_state_file"],
                    gate_cfg["sweep"]["initial_state_file"])
                self.assertNotIn("bias_points", gate_cfg["sweep"])

            vd005_drain = manifest["stages"][1]
            vd1_drain = manifest["stages"][4]
            self.assertAlmostEqual(
                0.005,
                vd005_drain["vela_physical_step_control"]["initial_step_V"])
            self.assertAlmostEqual(
                0.1,
                vd1_drain["vela_physical_step_control"]["initial_step_V"])

    def test_chain_validator_rejects_non_predecessor_restart(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m3_") as temporary:
            root = Path(temporary)
            manifest = MODULE.materialize(
                self.base(root), CONTRACT, root / "run", [0.05])
            manifest["stages"][2]["initial_state_file"] = "unrelated.csv"
            with self.assertRaisesRegex(ValueError, "restart state is not predecessor"):
                MODULE.validate_manifest_chain(manifest)

    @staticmethod
    def write_accepted_artifacts(stage: dict) -> None:
        config = MODULE.read_json(Path(stage["config"]))
        output = Path(stage["output_csv"])
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["bias_V", "current_total_A_per_um", "converged"])
            writer.writeheader()
            writer.writerow({
                "bias_V": config["sweep"]["stop"],
                "current_total_A_per_um": 0.0,
                "converged": 1,
            })
        Path(stage["final_state_file"]).write_text(
            "node_id,psi,phin,phip,electrons_m3,holes_m3\n"
            "0,0,0,0,1,1\n", encoding="utf-8")

    def test_execute_starts_successors_only_after_accepted_state(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m3_") as temporary:
            root = Path(temporary)
            output = root / "run"
            manifest = MODULE.materialize(
                self.base(root), CONTRACT, output, [0.05])
            calls: list[str] = []

            def fake_run(command, **_kwargs):
                config_path = Path(command[-1])
                stage = next(item for item in manifest["stages"]
                             if Path(item["config"]) == config_path)
                if stage["predecessor"] is not None:
                    predecessor = next(item for item in manifest["stages"]
                                       if item["id"] == stage["predecessor"])
                    self.assertEqual("accepted", predecessor["status"])
                    self.assertTrue(Path(stage["initial_state_file"]).is_file())
                    self.assertEqual(
                        predecessor["final_state_sha256"],
                        MODULE.sha256(Path(stage["initial_state_file"])))
                calls.append(stage["id"])
                self.write_accepted_artifacts(stage)
                return mock.Mock(returncode=0, stdout="ok\n", stderr="")

            with mock.patch.object(MODULE.subprocess, "run", side_effect=fake_run):
                result = MODULE.execute(manifest, root / "runner.exe", output)

            self.assertEqual("accepted", result["status"])
            self.assertEqual(3, len(calls))
            self.assertTrue(all(item["status"] == "accepted"
                                for item in result["stages"]))

    def test_rejected_stage_blocks_downstream_execution(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m3_") as temporary:
            root = Path(temporary)
            output = root / "run"
            manifest = MODULE.materialize(
                self.base(root), CONTRACT, output, [0.05])
            calls = 0

            def fake_run(_command, **_kwargs):
                nonlocal calls
                calls += 1
                stage = manifest["stages"][0]
                config = MODULE.read_json(Path(stage["config"]))
                Path(stage["output_csv"]).write_text(
                    "bias_V,converged\n"
                    f"{config['sweep']['stop']},0\n", encoding="utf-8")
                Path(stage["final_state_file"]).write_text("stale\n", encoding="utf-8")
                return mock.Mock(returncode=0, stdout="", stderr="")

            with mock.patch.object(MODULE.subprocess, "run", side_effect=fake_run):
                result = MODULE.execute(manifest, root / "runner.exe", output)

            self.assertEqual("fail", result["status"])
            self.assertEqual(1, calls)
            self.assertEqual("rejected", result["stages"][0]["status"])
            self.assertEqual("materialized", result["stages"][1]["status"])


if __name__ == "__main__":
    unittest.main()
