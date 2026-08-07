#!/usr/bin/env python3
"""Regression coverage for Sentaurus VM machine-readable diff report."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "generate_sentaurus_vm_diff_report.py"


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def write_common_exports(base: Path) -> None:
    write_csv(base / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0], [1, 1.0, 0.0]])
    write_csv(base / "elements.csv", ["id", "node0", "node1", "node2", "region_id"], [[0, 0, 1, 0, 0]])
    write_csv(base / "contacts.csv", ["id", "name", "region_id", "node_ids"], [[0, "anode", 0, "0;1"]])
    write_csv(base / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [[0, 1e17, 0.0], [1, 1e17, 0.0]])


def run_report(sentaurus: Path, vela: Path, tolerance: float = 1e-8, x_tolerance: float = 0.0) -> tuple[int, dict]:
    output = vela / "report.json"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--sentaurus-dir",
        str(sentaurus),
        "--vela-dir",
        str(vela),
        "--output-json",
        str(output),
        "--tolerance",
        str(tolerance),
        "--x-tolerance",
        str(x_tolerance),
    ]
    completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    payload = json.loads(output.read_text(encoding="utf-8"))
    return completed.returncode, payload


class GenerateSentaurusVmDiffReportTest(unittest.TestCase):
    def test_missing_candidate_and_missing_oracle_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_missing_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            write_common_exports(sentaurus)
            (sentaurus / "pn2d_sentaurus_inventory.json").write_text(json.dumps({"a": 1}), encoding="utf-8")
            write_csv(sentaurus / "reference_curves" / "pn2d_idvd_reference.csv", ["voltage_V", "current_A_per_um"], [[0.0, 0.0]])
            code, payload = run_report(sentaurus, vela)
            self.assertNotEqual(code, 0)
            self.assertFalse(payload["overall_pass"])
            reasons = {f["reason"] for f in payload["failures"]}
            self.assertIn("missing_candidate", reasons)

            sentaurus2 = root / "oracle2"
            vela2 = root / "vela2"
            sentaurus2.mkdir()
            vela2.mkdir()
            write_common_exports(vela2)
            (vela2 / "inventory.json").write_text(json.dumps({"a": 1}), encoding="utf-8")
            write_csv(vela2 / "pn2d_idvd.csv", ["voltage_V", "current_A_per_um"], [[0.0, 0.0]])
            code2, payload2 = run_report(sentaurus2, vela2)
            self.assertNotEqual(code2, 0)
            self.assertFalse(payload2["overall_pass"])
            reasons2 = {f["reason"] for f in payload2["failures"]}
            self.assertIn("missing_oracle", reasons2)

    def test_missing_both_curve_fails(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_missing_both_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            write_common_exports(sentaurus)
            write_common_exports(vela)
            inv = {"geometry": {"n": 1}}
            (sentaurus / "pn2d_sentaurus_inventory.json").write_text(json.dumps(inv), encoding="utf-8")
            (vela / "inventory.json").write_text(json.dumps(inv), encoding="utf-8")
            code, payload = run_report(sentaurus, vela)
            self.assertNotEqual(code, 0)
            self.assertIn({"item": "pn2d_idvd", "reason": "missing_both"}, payload["failures"])

    def test_curve_tolerance_and_x_tolerance(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_curve_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            write_common_exports(sentaurus)
            write_common_exports(vela)
            inv = {"geometry": {"n": 1}}
            (sentaurus / "pn2d_sentaurus_inventory.json").write_text(json.dumps(inv), encoding="utf-8")
            (vela / "inventory.json").write_text(json.dumps(inv), encoding="utf-8")
            write_csv(sentaurus / "reference_curves" / "pn2d_idvd_reference.csv", ["voltage_V", "current_A_per_um"], [[0.0, 0.0], [1.0, 1.0e-6]])
            write_csv(vela / "pn2d_idvd.csv", ["voltage_V", "current_A_per_um"], [[0.0, 0.0], [1.01, 1.08e-6]])

            code_pass, payload_pass = run_report(sentaurus, vela, tolerance=1e-7, x_tolerance=0.02)
            self.assertEqual(code_pass, 0)
            self.assertTrue(payload_pass["overall_pass"])

            code_y_fail, payload_y_fail = run_report(sentaurus, vela, tolerance=1e-8, x_tolerance=0.02)
            self.assertNotEqual(code_y_fail, 0)
            self.assertIn({"item": "pn2d_idvd", "reason": "curve_diff"}, payload_y_fail["failures"])

            code_x_fail, payload_x_fail = run_report(sentaurus, vela, tolerance=1e-7, x_tolerance=0.001)
            self.assertNotEqual(code_x_fail, 0)
            self.assertIn({"item": "pn2d_idvd", "reason": "curve_diff"}, payload_x_fail["failures"])

    def test_inventory_diff_and_invalid_curve_input_fail(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_inventory_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            write_common_exports(sentaurus)
            write_common_exports(vela)
            (sentaurus / "pn2d_sentaurus_inventory.json").write_text(json.dumps({"geometry": {"n": 1}}), encoding="utf-8")
            (vela / "inventory.json").write_text(json.dumps({"geometry": {"n": 2}}), encoding="utf-8")
            write_csv(sentaurus / "reference_curves" / "pn2d_idvd_reference.csv", ["voltage_V", "current_A_per_um"], [[0.0, 0.0], [1.0, 1.0e-6]])
            write_csv(vela / "pn2d_idvd.csv", ["voltage_V", "current_A_per_um"], [[0.0, 0.0], [1.0, 1.0e-6]])
            code, payload = run_report(sentaurus, vela)
            self.assertNotEqual(code, 0)
            self.assertIn({"item": "inventory", "reason": "json_diff"}, payload["failures"])

            write_csv(vela / "pn2d_idvd.csv", ["voltage_V", "current_A_per_um"], [[0.0, "nan"], [1.0, 1.0e-6]])
            output = vela / "report.json"
            cmd = [
                sys.executable,
                str(SCRIPT),
                "--sentaurus-dir",
                str(sentaurus),
                "--vela-dir",
                str(vela),
                "--output-json",
                str(output),
                "--tolerance",
                "1e-8",
                "--x-tolerance",
                "0.0",
            ]
            completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertFalse(payload["overall_pass"])
            self.assertEqual(payload["failures"][0]["reason"], "input_error")

    def test_decreasing_curve_passes_and_direction_mismatch_fails_with_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_direction_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            write_common_exports(sentaurus)
            write_common_exports(vela)
            inv = {"geometry": {"n": 1}}
            (sentaurus / "pn2d_sentaurus_inventory.json").write_text(json.dumps(inv), encoding="utf-8")
            (vela / "inventory.json").write_text(json.dumps(inv), encoding="utf-8")
            write_csv(sentaurus / "reference_curves" / "pn2d_idvd_reference.csv", ["voltage_V", "current_A_per_um"], [[2.0, 2.0e-6], [1.0, 1.0e-6], [0.0, 0.0]])
            write_csv(vela / "pn2d_idvd.csv", ["voltage_V", "current_A_per_um"], [[2.0, 2.0e-6], [1.0, 1.0e-6], [0.0, 0.0]])
            code_pass, payload_pass = run_report(sentaurus, vela)
            self.assertEqual(code_pass, 0)
            self.assertTrue(payload_pass["overall_pass"])

            write_csv(vela / "pn2d_idvd.csv", ["voltage_V", "current_A_per_um"], [[0.0, 0.0], [1.0, 1.0e-6], [2.0, 2.0e-6]])
            code_fail, payload_fail = run_report(sentaurus, vela)
            self.assertNotEqual(code_fail, 0)
            self.assertFalse(payload_fail["overall_pass"])
            self.assertEqual(payload_fail["failures"][0]["reason"], "input_error")

    def test_invalid_inventory_json_still_writes_failure_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_bad_json_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            write_common_exports(sentaurus)
            write_common_exports(vela)
            (sentaurus / "pn2d_sentaurus_inventory.json").write_text("{bad json", encoding="utf-8")
            (vela / "inventory.json").write_text(json.dumps({"geometry": {"n": 1}}), encoding="utf-8")
            code, payload = run_report(sentaurus, vela)
            self.assertNotEqual(code, 0)
            self.assertFalse(payload["overall_pass"])
            self.assertEqual(payload["failures"][0]["reason"], "input_error")

    def test_invalid_tolerance_precedes_missing_curve_and_replaces_old_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_bad_tolerance_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            output = vela / "report.json"
            output.write_text('{"stale": true}\n', encoding="utf-8")

            code, payload = run_report(sentaurus, vela, tolerance=-1.0)
            self.assertNotEqual(code, 0)
            self.assertFalse(payload["overall_pass"])
            self.assertEqual(payload["failures"][0]["item"], "arguments")
            self.assertEqual(payload["failures"][0]["reason"], "input_error")
            self.assertNotIn("stale", payload)

    def test_nonfinite_x_tolerance_writes_argument_failure_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_diff_bad_x_tolerance_") as tmp:
            root = Path(tmp)
            sentaurus = root / "oracle"
            vela = root / "vela"
            sentaurus.mkdir()
            vela.mkdir()
            for value in (float("nan"), float("inf")):
                code, payload = run_report(sentaurus, vela, x_tolerance=value)
                self.assertNotEqual(code, 0)
                self.assertFalse(payload["overall_pass"])
                self.assertEqual(payload["failures"][0]["item"], "arguments")
                self.assertEqual(payload["failures"][0]["reason"], "input_error")


if __name__ == "__main__":
    unittest.main()
