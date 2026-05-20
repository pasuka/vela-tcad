#!/usr/bin/env python3
"""Regression coverage for Sentaurus text import helpers."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class SentaurusImportToolsTest(unittest.TestCase):
    def test_plt_parser_exports_reference_curve(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_plt_") as tmp:
            root = Path(tmp)
            plt = root / "IdVd_n20_des.plt"
            out = root / "idvd_reference.csv"
            plt.write_text(
                """
DF-ISE text
Info {
  datasets = ["time" "gate OuterVoltage" "drain OuterVoltage" "drain TotalCurrent" "source TotalCurrent"]
}
Data {
  Values = [
    0 10 0 0 0
    1 10 30 2.19623e-3 -1.91281e-3
  ]
}
""".strip()
                + "\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "plt",
                    "--input",
                    str(plt),
                    "--output",
                    str(out),
                    "--bias-column",
                    "drain OuterVoltage",
                    "--current-column",
                    "drain TotalCurrent",
                ],
                check=True,
                cwd=REPO,
            )

            rows = out.read_text().splitlines()
            self.assertEqual(rows[0], "bias_V,current_total")
            self.assertIn("30,0.00219623", rows[2])

    def test_plt_parser_handles_multiline_data_block_without_values_label(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_plt_multiline_") as tmp:
            root = Path(tmp)
            plt = root / "IdVd_n20_des.plt"
            out = root / "idvd_reference.csv"
            plt.write_text(
                """
DF-ISE text
Info {
  datasets = ["time" "gate OuterVoltage" "drain OuterVoltage" "drain TotalCurrent"]
}
Data {
  0
  10 0
  0
  1
  10 30
  2.19623e-3
}
""".strip()
                + "\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "plt",
                    "--input",
                    str(plt),
                    "--output",
                    str(out),
                    "--bias-column",
                    "drain OuterVoltage",
                    "--current-column",
                    "drain TotalCurrent",
                ],
                check=True,
                cwd=REPO,
            )

            rows = out.read_text().splitlines()
            self.assertEqual(rows[0], "bias_V,current_total")
            self.assertEqual(rows[1], "0,0")
            self.assertEqual(rows[2], "30,0.00219623")

    def test_cmd_parser_exports_summary_and_python_runner(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_cmd_") as tmp:
            root = Path(tmp)
            cmd = root / "pp20_des.cmd"
            summary = root / "cmd_summary.json"
            runner = root / "run_idvd.py"
            cmd.write_text(
                """
File {
  grid = "n13_final_fps.tdr"
  current = "n20_des.plt"
  plot = "n20_des.tdr"
}
Electrode {
  { Name="drain" Voltage=0 }
  { Name="gate" Voltage=0 }
  { Name="source" Voltage=0 }
  { Name="substrate" Voltage=0 }
}
Physics { AreaFactor=15.0 Fermi Thermodynamic }
Solve {
  Coupled { Poisson Electron Hole }
  Quasistationary(
    Goal { Name="gate" Voltage=10 }
  ){ Coupled { Poisson Electron Hole } }
  NewCurrentFile="IdVd_"
  Quasistationary(
    Goal { Name="drain" Voltage=30 }
  ){ Coupled { Poisson Electron Hole } }
}
""".strip()
                + "\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "cmd",
                    "--input",
                    str(cmd),
                    "--summary-json",
                    str(summary),
                    "--python-runner",
                    str(runner),
                    "--mesh-json",
                    "vela/mesh.json",
                    "--output-csv",
                    "outputs/ldmos_idvd.csv",
                ],
                check=True,
                cwd=REPO,
            )

            data = json.loads(summary.read_text())
            self.assertEqual(data["files"]["grid"], "n13_final_fps.tdr")
            self.assertEqual(data["files"]["current"], "n20_des.plt")
            self.assertEqual({entry["name"] for entry in data["electrodes"]}, {"drain", "gate", "source", "substrate"})
            self.assertEqual(data["sweeps"][0]["contact"], "gate")
            self.assertEqual(data["sweeps"][0]["stop"], 10.0)
            self.assertEqual(data["sweeps"][1]["contact"], "drain")
            self.assertEqual(data["sweeps"][1]["stop"], 30.0)
            self.assertIn("Thermodynamic", data["unsupported_physics"])

            text = runner.read_text()
            self.assertIn("from vela.curves import run_iv_curve", text)
            self.assertIn("\"contact\": \"drain\"", text)
            self.assertIn("\"current_contact\": \"drain\"", text)


if __name__ == "__main__":
    unittest.main()
