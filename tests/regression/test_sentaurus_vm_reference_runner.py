#!/usr/bin/env python3
"""Regression coverage for the opt-in Sentaurus VM reference runner."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
RUNNER = REPO / "scripts" / "run_sentaurus_vm_reference.py"


class SentaurusVmReferenceRunnerTest(unittest.TestCase):
    def test_dry_run_writes_manifest_without_ssh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_vm_dry_") as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            for name in [
                "pn2d_sde.cmd",
                "pn2d_bv_sdevice.cmd",
                "models.par",
            ]:
                (source / name).write_text(f"{name}\n")
            out = root / "runs"

            subprocess.run([
                sys.executable,
                str(RUNNER),
                "pn2d",
                "--ssh-target", "sentaurus",
                "--source-dir", str(source),
                "--local-output-dir", str(out),
                "--remote-root", "~/sentaurus_runs/vela_oracle",
                "--run-id", "pn2d_bv_vm_dry_run",
                "--stages", "bv",
                "--dry-run",
            ], check=True)

            manifest = json.loads(
                (out / "pn2d_bv_vm_dry_run" / "sentaurus_vm_run_manifest.json").read_text()
            )
            self.assertEqual(manifest["ssh_target"], "sentaurus")
            self.assertEqual(
                manifest["remote_source_dir"],
                "~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source",
            )
            self.assertEqual(manifest["stages"], ["bv"])
            self.assertEqual(manifest["commands"], [
                "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source && sde -e -l pn2d_sde.cmd",
                "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source && sdevice pn2d_bv_sdevice.cmd > run_pn2d_bv.out 2>&1",
            ])

    def test_missing_required_deck_fails_before_ssh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_vm_missing_") as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            (source / "pn2d_sde.cmd").write_text("mesh\n")
            completed = subprocess.run([
                sys.executable,
                str(RUNNER),
                "pn2d",
                "--source-dir", str(source),
                "--local-output-dir", str(root / "runs"),
                "--run-id", "missing_bv",
                "--stages", "bv",
                "--dry-run",
            ], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("missing required source file", completed.stderr)


if __name__ == "__main__":
    unittest.main()
