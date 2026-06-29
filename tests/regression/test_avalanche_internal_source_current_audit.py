#!/usr/bin/env python3
"""Regression coverage for the avalanche internal source current audit runner."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_avalanche_internal_source_current_audit.py"


class AvalancheInternalSourceCurrentAuditRunnerTest(unittest.TestCase):
    def test_skip_run_writes_a2_b105_audit_config(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_internal_source_audit_") as td:
            tmp = Path(td)
            out_dir = tmp / "diagnostics"
            base_config = tmp / "base_bv.json"
            base_config.write_text(
                json.dumps({
                    "mesh_file": "mesh.json",
                    "node_doping_file": "doping.csv",
                    "materials_file": "materials.json",
                    "solver": {
                        "impact_ionization": {
                            "model": "van_overstraeten",
                            "generation": "current_density",
                        }
                    },
                    "sweep": {
                        "mode": "bv_reverse",
                        "contact": "Anode",
                        "current_contact": "Anode",
                        "bias_points": [0, -5],
                    },
                }),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--base-config", str(base_config),
                    "--out-dir", str(out_dir),
                    "--bias-points", "0,-5",
                    "--skip-run",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            config_path = (
                out_dir
                / "avalanche_internal_source_current_audit_case"
                / "BV-A2-B1p05-internal-source-current-audit"
                / "BV-A2-B1p05-internal-source-current-audit.json"
            )
            self.assertTrue(config_path.exists())
            config = json.loads(config_path.read_text(encoding="utf-8"))
            impact = config["solver"]["impact_ionization"]
            self.assertEqual(impact["model"], "van_overstraeten")
            self.assertEqual(impact["parameter_set"], "default")
            self.assertEqual(impact["driving_force"], "quasi_fermi_gradient")
            self.assertEqual(impact["generation"], "current_density")
            self.assertEqual(impact["A_scale"], 2.0)
            self.assertEqual(impact["B_scale"], 1.05)
            audit = config["sweep"]["diagnostics"]["avalanche_internal_source_current_audit"]
            self.assertTrue(audit["enabled"])
            self.assertEqual(
                Path(audit["csv_file"]),
                out_dir / "avalanche_internal_source_current_audit.csv",
            )
            self.assertEqual(
                Path(audit["summary_file"]),
                out_dir / "avalanche_internal_source_current_audit_summary.md",
            )


if __name__ == "__main__":
    unittest.main()
