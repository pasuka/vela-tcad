#!/usr/bin/env python3
"""Tests for the licensed-VM Sentaurus 2022 SDE fixture runner."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts import run_sentaurus2022_sde_fixture_vm as runner


REPO = Path(__file__).resolve().parents[2]
RUNNER = REPO / "scripts" / "run_sentaurus2022_sde_fixture_vm.py"


class Sentaurus2022SdeFixtureVmTest(unittest.TestCase):
    def test_wrapper_enables_evaluated_log_before_embedded_owned_source(self) -> None:
        self.assertEqual(
            runner.evaluated_wrapper("owned_sde.cmd", "(sde:clear)\n"),
            "(set! evaluate-log-file #t)\n"
            "; Vela owned source follows: owned_sde.cmd\n"
            "(sde:clear)\n",
        )

    def test_dry_run_seals_source_and_records_remote_commands(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_fixture_vm_") as tmp:
            root = Path(tmp)
            source = root / "owned_sde.cmd"
            source.write_text("(sde:clear)\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--source", str(source),
                    "--fixture-name", "owned",
                    "--local-output-dir", str(root / "out"),
                    "--remote-root", "~/sentaurus_runs/test",
                    "--run-id", "dry",
                    "--dry-run",
                ],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
            )
            manifest_path = root / "out" / "owned" / "dry" / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertTrue(manifest["dry_run"])
            self.assertEqual(manifest["fixture"], "owned")
            self.assertEqual(
                manifest["remote_run_dir"],
                "~/sentaurus_runs/test/owned_dry",
            )
            self.assertIn("input/owned_sde.cmd", manifest["files"])
            self.assertIn(
                "(set! evaluate-log-file #t)",
                (manifest_path.parent / "input" / "vela_sde_evaluated_wrapper.scm")
                .read_text(encoding="utf-8"),
            )
            self.assertIn(
                "(sde:clear)",
                (manifest_path.parent / "input" / "vela_sde_evaluated_wrapper.scm")
                .read_text(encoding="utf-8"),
            )
            self.assertEqual(json.loads(completed.stdout)["source_sha256"], manifest["source_sha256"])

    def test_mesh_summary_uses_final_mesh_counts(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_mesh_summary_") as tmp:
            stdout = Path(tmp) / "sde.stdout.txt"
            stdout.write_text(
                "The final mesh has:\n    1943 points.\n"
                "    3680 elements.\n       3 regions.\n",
                encoding="utf-8",
            )
            self.assertEqual(
                runner.mesh_summary(stdout),
                {"points": 1943, "elements": 3680, "regions": 3},
            )


if __name__ == "__main__":
    unittest.main()
