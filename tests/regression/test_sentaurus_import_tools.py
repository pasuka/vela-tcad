#!/usr/bin/env python3
"""Regression coverage for Sentaurus text import helpers."""

from __future__ import annotations

import json
import os
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
            deck = root / "run_idvd_deck.json"
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
Physics (MaterialInterface="Oxide/Silicon") { Traps(Conc=5e8) }
Plot {
  TotalCurrent/Vector ElectricField/Vector Potential Doping
}
Math {
  Extrapolate
  Iterations=12
}
Solve {
  Poisson
  Coupled(Iterations=100){ Poisson Electron }
  Coupled { Poisson Electron Hole }
  Quasistationary(
    InitialStep=1e-4 MinStep=1e-5 MaxStep=1
    Goal { Name="gate" Voltage=10 }
  ){ Coupled { Poisson Electron Hole } }
  NewCurrentFile="IdVd_"
  Quasistationary(
    InitialStep=1e-6 MinStep=1e-12 MaxStep=0.1
    Goal { Name="drain" Voltage=30 }
  ){ Coupled { Poisson Electron Hole }
     CurrentPlot( Time=(Range=(0 1) Intervals=80 ) )
  }
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
                    "--deck-json",
                    str(deck),
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
            self.assertEqual(data["physics"][0]["scope"], {"kind": "global"})
            self.assertEqual(data["physics"][0]["parameters"]["AreaFactor"], 15.0)
            self.assertIn("Fermi", data["physics"][0]["models"])
            self.assertIn(
                {"name": "TotalCurrent", "vector": True},
                data["plot_fields"],
            )
            self.assertEqual(data["math"]["parameters"]["Iterations"], 12.0)
            self.assertIn("Extrapolate", data["math"]["flags"])
            self.assertEqual(data["solve"]["initial_steps"][1]["equations"], ["Poisson", "Electron"])
            self.assertEqual(data["sweeps"][1]["step_control"]["InitialStep"], 1.0e-6)
            self.assertEqual(data["sweeps"][1]["current_plot"]["intervals"], 80)
            self.assertIn("Thermodynamic", data["unsupported_physics"])
            self.assertIn("Traps", data["unsupported_physics"])
            self.assertIn("Thermodynamic", [item["feature"] for item in data["unsupported_report"]])
            self.assertIn("Traps", [item["feature"] for item in data["unsupported_report"]])

            text = runner.read_text()
            self.assertIn("from vela.curves import run_iv_curve", text)
            self.assertIn("\"contact\": \"drain\"", text)
            self.assertIn("\"current_contact\": \"drain\"", text)
            self.assertIn("\"unsupported_report\"", text)
            deck_data = json.loads(deck.read_text())
            self.assertEqual(deck_data["sweep"]["contact"], "drain")
            self.assertEqual(deck_data["sentaurus_import"]["sweeps"][0]["contact"], "gate")

    def test_cmd_parser_expands_template_variables(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_cmd_template_") as tmp:
            root = Path(tmp)
            cmd = root / "IdVd_des.cmd"
            summary = root / "cmd_summary.json"
            cmd.write_text(
                """
File {
  grid = "n@node|DevMesh@_final_fps.tdr"
  current = "@plot@"
  plot = "@tdrdat@"
}
Electrode {
  { Name="drain" Voltage=0.0 }
  { Name="gate" Voltage=0.0 }
}
Thermode {
  { Name="drain" Temperature=300 SurfaceResistance=@ThermalR@ }
}
Physics { AreaFactor=15.0 Thermodynamic }
Plot { ElectricField/Vector Potential Doping }
Math { Iterations=12 ExitOnFailure }
Solve {
  Quasistationary(
    InitialStep=1e-4 MinStep=1e-5 MaxStep=1
    Goal { Name="gate" Voltage=@Vg@ }
  ){ Coupled { Poisson Electron Hole Temperature } }
  NewCurrentFile="IdVd_"
  Quasistationary(
    InitialStep=1e-6 MinStep=1e-12 MaxStep=0.1
    Goal { Name="drain" Voltage=30.0 }
  ){ Coupled { Poisson Electron Hole Temperature }
     CurrentPlot( Time=(Range=(0 1) Intervals=80 ) )
  }
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
                    "--template-var",
                    "node|DevMesh=13",
                    "--template-var",
                    "plot=n20_des.plt",
                    "--template-var",
                    "tdrdat=n20_des.tdr",
                    "--template-var",
                    "ThermalR=5.0e-4",
                    "--template-var",
                    "Vg=10.0",
                ],
                check=True,
                cwd=REPO,
            )

            data = json.loads(summary.read_text())
            self.assertEqual(data["files"]["grid"], "n13_final_fps.tdr")
            self.assertEqual(data["files"]["current"], "n20_des.plt")
            self.assertEqual(data["files"]["plot"], "n20_des.tdr")
            self.assertEqual(data["template_variables"]["Vg"], "10.0")
            self.assertEqual(data["unresolved_placeholders"], [])
            self.assertEqual(data["thermodes"][0]["surface_resistance"], 5.0e-4)
            self.assertEqual(data["sweeps"][0]["contact"], "gate")
            self.assertEqual(data["sweeps"][0]["stop"], 10.0)
            self.assertEqual(data["sweeps"][0]["equations"], ["Poisson", "Electron", "Hole", "Temperature"])
            self.assertEqual(data["sweeps"][1]["current_plot"]["intervals"], 80)

    def test_project_import_generates_neutral_reference_tree(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_project_") as tmp:
            root = Path(tmp)
            project = root / "project"
            output = root / "sentaurus_ldmos_n20"
            project.mkdir()
            (project / "n20_des.tdr").write_text("fake tdr placeholder\n")
            (project / "IdVd_n20_des.plt").write_text(
                """
DF-ISE text
Info {
  datasets = ["time" "gate OuterVoltage" "drain OuterVoltage" "drain TotalCurrent"]
}
Data {
  0 10 0 0
  1 10 30 2.19623e-3
}
""".strip()
                + "\n"
            )
            (project / "pp20_des.cmd").write_text(
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
Physics { Thermodynamic }
Solve {
  Quasistationary( Goal { Name="gate" Voltage=10 } ){ Coupled { Poisson Electron Hole } }
  Quasistationary( Goal { Name="drain" Voltage=30 } ){ Coupled { Poisson Electron Hole } }
}
""".strip()
                + "\n"
            )
            fake_importer = root / "fake_tdr_importer.py"
            fake_importer.write_text(
                """
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--tdr", required=True)
parser.add_argument("--inventory-json")
parser.add_argument("--export-dir")
args = parser.parse_args()

if args.export_dir:
    out = Path(args.export_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "fields").mkdir(exist_ok=True)
    def write_csv(path, header, rows):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)
    write_csv(out / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0, 0], [1, 1, 0], [2, 1, 1], [3, 0, 1]])
    write_csv(out / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [[0, 0, 1, 2, "Silicon_1", "Si"], [1, 0, 2, 3, "Silicon_1", "Si"]])
    write_csv(out / "contacts.csv", ["name", "node_ids", "region"], [["drain", "1;2", "Silicon_1"], ["gate", "0;3", "Silicon_1"]])
    write_csv(out / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [[0, 1e17, 0], [1, 1e17, 0], [2, 1e17, 0], [3, 1e17, 0]])
    (out / "metadata.json").write_text(json.dumps({"vertex_count": 4, "region_count": 1, "dataset_count": 0}, indent=2) + "\\n")
    (out / "field_manifest.json").write_text(json.dumps({"fields": []}, indent=2) + "\\n")

if args.inventory_json:
    Path(args.inventory_json).write_text(json.dumps({"vertex_count": 4, "region_count": 1, "dataset_count": 0}, indent=2) + "\\n")
""".strip()
                + "\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "project",
                    "--project-dir",
                    str(project),
                    "--node",
                    "20",
                    "--device",
                    "ldmos2d",
                    "--output-dir",
                    str(output),
                    "--tdr-importer",
                    f"{sys.executable} {fake_importer}",
                ],
                check=True,
                cwd=REPO,
                env={**os.environ, "PYTHONIOENCODING": "utf-8"},
            )

            for required in [
                "nodes.csv",
                "elements.csv",
                "contacts.csv",
                "doping.csv",
                "field_manifest.json",
                "metadata.json",
                "import_summary.json",
                "cmd_summary.json",
                "run_idvd.py",
                "run_idvd_deck.json",
                "reference_curves/ldmos2d_idvd_reference.csv",
                "vela/mesh.json",
                "vela/simulation_iv.json",
            ]:
                self.assertTrue((output / required).is_file(), required)

            summary = json.loads((output / "import_summary.json").read_text())
            self.assertEqual(summary["node"], 20)
            self.assertEqual(summary["device"], "ldmos2d")
            self.assertEqual(summary["sources"]["tdr"], str(project / "n20_des.tdr"))
            self.assertIn("Thermodynamic", summary["warnings"][0])
            self.assertIn("field_manifest.json", summary["generated"])
            self.assertIn("reference_curves/ldmos2d_idvd_reference.csv", summary["generated"])


if __name__ == "__main__":
    unittest.main()
