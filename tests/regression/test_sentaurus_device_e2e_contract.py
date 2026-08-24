#!/usr/bin/env python3
"""Unified SDE/SDevice deck and checked-in acceptance-contract tests."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "sentaurus_import", REPO / "scripts" / "sentaurus_import.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

from sentaurus_sde_frontend import (  # noqa: E402
    offline_static_frontend,
    sentaurus2022_evaluated_frontend,
)


def execution(kind: str, sweep_count: int = 1) -> dict:
    stages = [{
        "index": 0, "phase": "initial", "type": "Coupled",
        "equations": ["Poisson"], "depends_on": [],
    }]
    for index in range(sweep_count):
        stages.append({
            "index": index + 1, "phase": "sweep", "type": "Quasistationary",
            "equations": ["Poisson", "Electron", "Hole"],
            "depends_on": [index], "contact": "Anode",
            "goal_voltage": -0.2 if kind == "bv" else 0.3,
            "step_control": {},
        })
    return {
        "source": "case_des.cmd",
        "electrodes": [{"name": "Anode", "voltage": 0.0},
                       {"name": "Cathode", "voltage": 0.0}],
        "stages": stages,
        "analysis": {"kind": kind, "final_contact": "Anode",
                     "final_voltage": -0.2 if kind == "bv" else 0.3},
    }


class UnifiedDeckTest(unittest.TestCase):
    DEVICE = {"source": "case_sde.cmd", "contacts": [
        {"name": "Anode"}, {"name": "Cathode"}]}

    def test_device_cli_defaults_to_license_free_offline_frontend(self) -> None:
        args = MODULE.build_parser().parse_args([
            "device",
            "--sde", "case_sde.cmd",
            "--sdevice", "case_sdevice.cmd",
            "--output-dir", "generated",
        ])
        self.assertEqual("offline-static", args.sde_frontend)

    def test_bv_mode_and_ordered_execution_sequence(self) -> None:
        deck, _warnings = MODULE.build_device_deck(
            self.DEVICE, execution("bv"), {"physics": []},
            "mesh.json", "doping.csv", "bv.csv")
        self.assertEqual("bv_reverse", deck["sweep"]["mode"])
        self.assertEqual([0, 1], [stage["index"]
                                  for stage in deck["_sentaurus_execution_sequence"]])
        self.assertEqual("gummel_on_newton_failure",
                         deck["solver"]["handoff"]["fallback"])
        self.assertFalse(
            deck["solver"]["handoff"]["require_gummel_convergence"])
        self.assertEqual(
            "report", deck["solver"]["global_continuity_closure"]["mode"])
        self.assertEqual(
            0.05, deck["solver"]["global_continuity_closure"]["tolerance"])
        self.assertLess(deck["sweep"]["shrink_factor"], 1.0)

    def test_quasistationary_step_controls_scale_to_voltage_span(self) -> None:
        ir = execution("iv")
        ir["stages"][-1]["step_control"] = {
            "InitialStep": 1.0e-3,
            "MinStep": 1.0e-8,
            "MaxStep": 0.025,
            "Increment": 1.3,
            "Decrement": 2.0,
        }
        deck, _warnings = MODULE.build_device_deck(
            self.DEVICE, ir, {"physics": []},
            "mesh.json", "doping.csv", "iv.csv")
        self.assertAlmostEqual(3.0e-4, deck["sweep"]["initial_step"])
        self.assertAlmostEqual(3.0e-9, deck["sweep"]["min_step"])
        self.assertAlmostEqual(7.5e-3, deck["sweep"]["max_step"])
        self.assertAlmostEqual(1.3, deck["sweep"]["growth_factor"])
        self.assertAlmostEqual(0.5, deck["sweep"]["shrink_factor"])
        self.assertEqual(
            "gummel_on_newton_failure", deck["solver"]["handoff"]["fallback"])
        self.assertEqual(300, deck["solver"]["handoff"]["gummel_max_iter"])
        self.assertTrue(
            deck["solver"]["handoff"]["require_gummel_convergence"])
        self.assertEqual(
            "componentwise", deck["solver"]["gummel_carrier_convergence_norm"])
        self.assertEqual(
            "report", deck["solver"]["carrier_row_convergence"]["mode"])
        self.assertEqual(
            "report", deck["solver"]["global_continuity_closure"]["mode"])

    def test_multiple_sweeps_fail_closed_instead_of_dropping_earlier_sweep(self) -> None:
        with self.assertRaises(MODULE.ExecutionIrError) as ctx:
            MODULE.build_device_deck(
                self.DEVICE, execution("iv", 2), {"physics": []},
                "mesh.json", "doping.csv", "iv.csv")
        self.assertIn("refusing to silently discard stages", str(ctx.exception))

    def test_qualified_bv_profile_restores_validated_newton_policy(self) -> None:
        deck, _warnings = MODULE.build_device_deck(
            self.DEVICE, execution("bv"), {"physics": []},
            "mesh.json", "doping.csv", "bv.csv",
            solver_profile="pn2d_sentaurus2018_bv_qualified")
        solver = deck["solver"]
        self.assertEqual("pn2d_sentaurus2018_bv_qualified",
                         deck["_e2e_solver_profile"])
        self.assertEqual(0, solver["handoff"]["gummel_max_iter"])
        self.assertEqual("none", solver["handoff"]["fallback"])
        self.assertFalse(solver["handoff"]["require_gummel_convergence"])
        self.assertFalse(solver["mobility"]["jacobian_field_derivatives"])
        self.assertNotIn("quasi_fermi_reference", solver)
        self.assertEqual("report", solver["global_continuity_closure"]["mode"])
        self.assertEqual(0.05, deck["sweep"]["initial_step"])
        self.assertEqual(0.05, deck["sweep"]["max_step"])

    def test_unknown_solver_profile_fails_closed(self) -> None:
        with self.assertRaises(MODULE.ExecutionIrError):
            MODULE.build_device_deck(
                self.DEVICE, execution("bv"), {"physics": []},
                "mesh.json", "doping.csv", "bv.csv",
                solver_profile="not-a-profile")


class ContractGateTest(unittest.TestCase):
    def test_full_bv_contract_targets_original_minus_twenty_volts(self) -> None:
        contract = MODULE.load_e2e_contract()
        case = contract["cases"]["pn2d_sentaurus2018_bv"]
        self.assertEqual("bv", case["analysis"])
        self.assertEqual("pn2d_sentaurus2018_bv_qualified",
                         case["solver_profile"])
        self.assertEqual(-20.0, case["curve"]["bias_stop"])
        self.assertGreaterEqual(case["curve"]["min_points"], 20)

    def test_partial_curve_reports_last_converged_bias(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_partial_curve_") as tmp:
            curve = Path(tmp) / "bv.csv"
            curve.write_text(
                "bias_V,converged\n0,1\n-0.5,1\n-1.0,0\n",
                encoding="utf-8")
            self.assertEqual(-0.5, MODULE._last_converged_bias(curve))

    def test_checked_in_contract_is_consumed_and_can_accept_a_curve(self) -> None:
        contract = MODULE.load_e2e_contract()
        case = contract["cases"]["pn2d_sentaurus2018_coarse7x3_iv"]
        qualification = {
            "node_count": 21, "cell_count": 24, "max_angle_degrees": 90.0,
            "node_volume_policy": "mixed_voronoi",
            "contacts": [{"node_count": 3}, {"node_count": 3}],
        }
        with tempfile.TemporaryDirectory(prefix="vela_contract_") as tmp:
            curve = Path(tmp) / "iv.csv"
            with curve.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "mode", "bias_V", "converged", "current_total_A_per_um",
                    "global_electron_continuity_closure_ratio",
                    "global_hole_continuity_closure_ratio",
                ])
                writer.writeheader()
                for index in range(10):
                    writer.writerow({
                        "mode": "iv", "bias_V": 0.3 * index / 9.0,
                        "converged": 1, "current_total_A_per_um": index * 1e-9,
                        "global_electron_continuity_closure_ratio": 0.01,
                        "global_hole_continuity_closure_ratio": 0.01,
                    })
            terminal = Path(tmp) / "terminal_currents.csv"
            with terminal.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "bias_V", "contact", "I_sgflux_A_per_um"])
                writer.writeheader()
                for index in range(10):
                    bias = 0.3 * index / 9.0
                    writer.writerow({"bias_V": bias, "contact": "Anode",
                                     "I_sgflux_A_per_um": index * 1e-9})
                    writer.writerow({"bias_V": bias, "contact": "Cathode",
                                     "I_sgflux_A_per_um": -index * 1e-9})
            report = MODULE.evaluate_e2e_acceptance(
                "pn2d_sentaurus2018_coarse7x3_iv", case, "iv",
                qualification, curve, terminal)
        self.assertEqual("accepted", report["status"])

    def test_mixed_terminal_current_tolerance_handles_near_zero_current(self) -> None:
        contract = MODULE.load_e2e_contract()
        case = contract["cases"]["pn2d_sentaurus2018_iv"]
        qualification = {
            "node_count": 1767, "cell_count": 3360,
            "max_angle_degrees": 90.0,
            "node_volume_policy": "mixed_voronoi",
            "contacts": [{"node_count": 31}, {"node_count": 31}],
        }
        with tempfile.TemporaryDirectory(prefix="vela_mixed_current_gate_") as tmp:
            curve = Path(tmp) / "iv.csv"
            terminal = Path(tmp) / "terminal.csv"
            with curve.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "mode", "bias_V", "converged", "current_total_A_per_um"])
                writer.writeheader()
                for index in range(10):
                    writer.writerow({
                        "mode": "iv", "bias_V": 0.3 * index / 9,
                        "converged": 1,
                        "current_total_A_per_um": index * 1.0e-16,
                    })
            with terminal.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "bias_V", "contact", "I_sgflux_A_per_um"])
                writer.writeheader()
                for index in range(10):
                    bias = 0.3 * index / 9
                    current = max(index, 1) * 1.0e-20
                    writer.writerow({"bias_V": bias, "contact": "Anode",
                                     "I_sgflux_A_per_um": current})
                    writer.writerow({"bias_V": bias, "contact": "Cathode",
                                     "I_sgflux_A_per_um": -0.5 * current})
            report = MODULE.evaluate_e2e_acceptance(
                "pn2d_sentaurus2018_iv", case, "iv", qualification,
                curve, terminal)
        self.assertEqual("accepted", report["status"])
        self.assertGreater(report["checks"]["max_terminal_current_imbalance"], 0.02)
        self.assertLess(
            report["checks"]["max_terminal_current_mixed_tolerance_ratio"], 1.0)

    def test_matching_contract_clamps_generated_sweep_to_acceptance_stop(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_contract_deck_") as tmp:
            subprocess.run([
                sys.executable, str(REPO / "scripts" / "sentaurus_import.py"),
                "device", "--sde",
                str(REPO / "reference_tcad" / "pn2d_sentaurus2018_coarse7x3" /
                    "source" / "pn2d_sde.cmd"),
                "--sdevice",
                str(REPO / "reference_tcad" / "pn2d_sentaurus2018" /
                    "source" / "pn2d_iv_sdevice.cmd"),
                "--output-dir", tmp,
            ], cwd=REPO, check=True, capture_output=True, text=True)
            deck = json.loads((Path(tmp) / "simulation.json").read_text())
            manifest = json.loads((Path(tmp) / "run_manifest.json").read_text())
        self.assertEqual(0.3, deck["sweep"]["stop"])
        self.assertEqual(0.01, deck["sweep"]["initial_step"])
        self.assertEqual(0.25, deck["sweep"]["max_step"])
        self.assertEqual("e2e_contract", manifest["sweep_stop_source"])


class SdeFrontendProvenanceTest(unittest.TestCase):
    def test_offline_manifest_has_source_hash(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_frontend_") as tmp:
            source = Path(tmp) / "case_sde.cmd"
            source.write_text(
                '(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) '
                '"Silicon" "R.Si")\n', encoding="utf-8")
            _ir, manifest = offline_static_frontend(source)
        self.assertEqual("offline-static", manifest["frontend"])
        self.assertEqual(64, len(manifest["source_sha256"]))
        self.assertIsNone(manifest["official_mesh_tdr"])

    def test_evaluated_frontend_builds_wrapper_and_hashes_log(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_evaluated_") as tmp:
            root = Path(tmp)
            source = root / "case_sde.cmd"
            source.write_text('(define L 1)\n', encoding="utf-8")
            output = root / "out"
            calls = 0

            def fake_run(command, **_kwargs):
                nonlocal calls
                calls += 1
                if calls == 1:
                    self.assertEqual(command, ["sde", "-v"])
                    return SimpleNamespace(stdout="SDE T-2022.03", stderr="", returncode=0)
                (output / "vela_sde_evaluated_wrapper.log.eval").write_text(
                    '(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) '
                    '"Silicon" "R.Si")\n', encoding="utf-8")
                return SimpleNamespace(stdout="ok", stderr="", returncode=0)

            with patch("sentaurus_sde_frontend.subprocess.run", side_effect=fake_run):
                ir, manifest = sentaurus2022_evaluated_frontend(
                    source, output, executable="sde")
            wrapper = (output / "vela_sde_evaluated_wrapper.scm").read_text()
        self.assertIn("(set! evaluate-log-file #t)", wrapper)
        self.assertIn("(define L 1)", wrapper)
        self.assertNotIn('(load "', wrapper)
        self.assertEqual("sentaurus2022-evaluated", manifest["frontend"])
        self.assertEqual("SDE T-2022.03", manifest["sde_version"])
        self.assertEqual(64, len(manifest["evaluated_log_sha256"]))
        self.assertIsNone(manifest["official_mesh_tdr"])
        self.assertEqual("R.Si", ir["geometry"]["regions"][0]["name"])

    def test_evaluated_frontend_records_official_mesh_tdr(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_evaluated_mesh_") as tmp:
            root = Path(tmp)
            source = root / "case_sde.cmd"
            source.write_text('(define L 1)\n', encoding="utf-8")
            output = root / "out"
            calls = 0

            def fake_run(_command, **_kwargs):
                nonlocal calls
                calls += 1
                if calls == 1:
                    return SimpleNamespace(
                        stdout="SDE T-2022.03", stderr="", returncode=0)
                (output / "case.log.eval").write_text(
                    '(sdegeo:create-rectangle (position 0 0 0) '
                    '(position 1 1 0) "Silicon" "R.Si")\n',
                    encoding="utf-8")
                (output / "sentaurus_official_msh.tdr").write_bytes(b"tdr")
                return SimpleNamespace(stdout="ok", stderr="", returncode=0)

            with patch("sentaurus_sde_frontend.subprocess.run", side_effect=fake_run):
                _ir, manifest = sentaurus2022_evaluated_frontend(
                    source, output, executable="sde", build_mesh=True)
        self.assertTrue(str(manifest["official_mesh_tdr"]).endswith("_msh.tdr"))
        self.assertEqual(64, len(manifest["official_mesh_tdr_sha256"]))


if __name__ == "__main__":
    unittest.main()
