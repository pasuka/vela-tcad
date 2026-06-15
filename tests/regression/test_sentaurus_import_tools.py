#!/usr/bin/env python3
"""Regression coverage for Sentaurus text import helpers."""

from __future__ import annotations

import importlib.util
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SENTAURUS_IMPORT_SPEC = importlib.util.spec_from_file_location(
    "sentaurus_import",
    REPO / "scripts" / "sentaurus_import.py",
)
assert SENTAURUS_IMPORT_SPEC is not None
sentaurus_import = importlib.util.module_from_spec(SENTAURUS_IMPORT_SPEC)
assert SENTAURUS_IMPORT_SPEC.loader is not None
SENTAURUS_IMPORT_SPEC.loader.exec_module(sentaurus_import)


class SentaurusImportToolsTest(unittest.TestCase):
    def test_solver_physics_ignores_bare_mobility_model(self) -> None:
        deck = {"solver": {"type": "gummel"}}
        warnings = sentaurus_import.apply_solver_physics(
            deck,
            {"physics": [{"models": ["Mobility"]}]},
            {"name": "iv", "kind": "iv"},
        )

        self.assertNotIn("mobility", deck["solver"])
        self.assertEqual(warnings, [])

    def test_solver_physics_maps_doping_dependence_to_masetti(self) -> None:
        deck = {"solver": {"type": "gummel"}}

        warnings = sentaurus_import.apply_solver_physics(
            deck,
            {"physics": [{"models": ["Mobility", "DopingDependence"]}]},
            {"name": "iv", "kind": "iv"},
        )

        self.assertEqual(deck["solver"]["mobility"]["model"], "masetti")
        self.assertEqual(warnings, [])

    def test_solver_physics_maps_doping_dependence_high_field_to_masetti_field(self) -> None:
        deck = {"solver": {"type": "gummel"}}

        warnings = sentaurus_import.apply_solver_physics(
            deck,
            {"physics": [{"models": ["Mobility", "DopingDependence", "HighFieldSaturation"]}]},
            {"name": "iv", "kind": "iv"},
        )

        self.assertEqual(deck["solver"]["mobility"]["model"], "masetti_field")
        self.assertEqual(warnings, [])

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

    def test_sde_parser_exports_pn2d_summary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_sde_") as tmp:
            root = Path(tmp)
            sde = root / "pn2d_sde.cmd"
            summary = root / "sde_summary.json"
            sde.write_text(
                """
(define L  2.0)
(define H  0.5)
(define Xj 1.0)
(define Na 1e17)
(define Nd 1e17)
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position L H 0.0)
  "Silicon"
  "R.Si")
(sdedr:define-constant-profile "P_Doping" "BoronActiveConcentration" Na)
(sdedr:define-constant-profile-region "P_Doping_Reg" "P_Doping" "R.Si")
(sdedr:define-constant-profile "N_Doping" "PhosphorusActiveConcentration" Nd)
(sdegeo:create-rectangle
  (position Xj 0.0 0.0)
  (position L H 0.0)
  "Silicon"
  "R.NRegion")
(sdedr:define-constant-profile-region "N_Doping_Reg" "N_Doping" "R.NRegion")
(sdegeo:define-contact-set "Anode" 4.0 (color:rgb 1 0 0) "##")
(sdegeo:define-contact-set "Cathode" 4.0 (color:rgb 0 0 1) "##")
(sdedr:define-refinement-window "JunctionWindow"
  "Rectangle"
  (position (- Xj 0.15) 0.0 0.0)
  (position (+ Xj 0.15) H 0.0))
""".strip()
                + "\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "sde",
                    "--input",
                    str(sde),
                    "--summary-json",
                    str(summary),
                ],
                check=True,
                cwd=REPO,
            )

            data = json.loads(summary.read_text())
            self.assertEqual(data["defines"]["L"], 2.0)
            self.assertEqual(data["defines"]["H"], 0.5)
            self.assertEqual(data["defines"]["Xj"], 1.0)
            self.assertEqual(data["geometry"]["rectangles"][0]["region"], "R.Si")
            self.assertEqual(data["doping_profiles"]["P_Doping"]["species"], "BoronActiveConcentration")
            self.assertEqual(data["doping_profiles"]["P_Doping"]["value"], 1.0e17)
            self.assertEqual(data["doping_placements"]["P_Doping_Reg"]["region"], "R.Si")
            self.assertEqual({item["name"] for item in data["contacts"]}, {"Anode", "Cathode"})
            self.assertEqual(data["refinement_windows"][0]["name"], "JunctionWindow")
            self.assertEqual(data["refinement_windows"][0]["lower_left"], [0.85, 0.0])
            self.assertEqual(data["refinement_windows"][0]["upper_right"], [1.15, 0.5])

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
                "reference_tcad_manifest.json",
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
            self.assertTrue(Path(summary["sources"]["tdr"]).samefile(project / "n20_des.tdr"))
            self.assertIn("Thermodynamic", summary["warnings"][0])
            self.assertIn("field_manifest.json", summary["generated"])
            self.assertIn("reference_tcad_manifest.json", summary["generated"])
            self.assertIn("reference_curves/ldmos2d_idvd_reference.csv", summary["generated"])

            manifest = json.loads((output / "reference_tcad_manifest.json").read_text())
            self.assertEqual(manifest["schema"], "vela.reference_tcad.sentaurus_project.v1")
            self.assertEqual(manifest["device"], "ldmos2d")
            self.assertEqual(manifest["node"], 20)
            self.assertFalse(manifest["commit_policy"]["raw_sentaurus_artifacts"])
            self.assertIn("reference_curves/ldmos2d_idvd_reference.csv", manifest["reference_curves"])
            self.assertIn("vela/simulation_iv.json", manifest["vela_decks"])
            self.assertIn("unsupported physics: Thermodynamic", manifest["warnings"])

    def test_reference_import_config_supports_equilibrium_newton_fixture(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_equilibrium_reference_") as tmp:
            root = Path(tmp)
            source = root / "pn2d_sentaurus2018"
            output = root / "out"
            source.mkdir()
            for name in ["pn2d_msh.tdr", "pn2d_0v_des.tdr", "pn2d_iv_des.tdr"]:
                (source / name).write_text("fake tdr placeholder\n")
            (source / "pn2d_sde.cmd").write_text("(define L 2.0)\n(define H 0.5)\n(define XJ 1.0)\n")
            (source / "pn2d_0v_sdevice.cmd").write_text(
                """
File { Grid="pn2d_msh.tdr" Plot="pn2d_0v_des.tdr" Current="pn2d_0v.plt" }
Electrode { { Name="Anode" Voltage=0 } { Name="Cathode" Voltage=0 } }
Physics { Mobility(DopingDependence) Recombination(SRH) EffectiveIntrinsicDensity(OldSlotboom) }
Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled(Iterations=100) { Poisson Electron Hole }
}
""".strip()
                + "\n"
            )
            (source / "pn2d_iv_sdevice.cmd").write_text(
                """
File { Grid="pn2d_msh.tdr" Plot="pn2d_iv_des.tdr" Current="pn2d_iv.plt" }
Electrode { { Name="Anode" Voltage=0 } { Name="Cathode" Voltage=0 } }
Physics { Mobility(DopingDependence) Recombination(SRH) EffectiveIntrinsicDensity(OldSlotboom) }
Solve {
  Coupled(Iterations=100) { Poisson Electron Hole }
  Quasistationary( MaxStep=0.05 Goal { Name="Anode" Voltage=1.0 } ){ Coupled { Poisson Electron Hole } }
}
""".strip()
                + "\n"
            )
            (source / "pn2d_0v.plt").write_text(
                """
Info { datasets = ["time" "Anode OuterVoltage" "Anode TotalCurrent" "Cathode TotalCurrent"] }
Data {
  0 0 -2.0e-25 2.0e-25
  0 0 -7.0e-25 7.0e-25
}
""".strip()
                + "\n"
            )
            (source / "pn2d_iv.plt").write_text(
                """
Info { datasets = ["time" "Anode OuterVoltage" "Anode TotalCurrent" "Cathode TotalCurrent"] }
Data {
  0 0 0 0
  1 1.0 1.25e-4 -1.25e-4
}
""".strip()
                + "\n"
            )
            config = root / "pn2d_sentaurus2018_reference.json"
            config.write_text(json.dumps({
                "case": "pn2d_sentaurus2018",
                "device": "pn_diode",
                "mesh_tdr": "pn2d_msh.tdr",
                "sde_cmd": "pn2d_sde.cmd",
                "tdr_doping": {"compensated_node_policy": "dominant_signed_region"},
                "vela_solver": {
                    "method": "gummel_newton",
                    "max_iter": 40,
                    "reltol": 1.0e-8,
                    "abstol": 1.0e-18,
                    "max_update": 5.0,
                    "line_search": True,
                    "warm_start": True,
                    "handoff": {
                        "fallback": "none",
                        "require_gummel_convergence": False,
                        "gummel_max_iter": 0,
                        "newton_max_iter": 40,
                    },
                },
                "simulations": [
                    {
                        "name": "0v",
                        "kind": "equilibrium",
                        "tdr": "pn2d_0v_des.tdr",
                        "cmd": "pn2d_0v_sdevice.cmd",
                        "plt": "pn2d_0v.plt",
                        "bias_column": "Anode OuterVoltage",
                        "current_column": "Anode TotalCurrent",
                        "vela_solver": {
                            "max_iter": 80,
                            "reltol": 1.0e-10,
                            "abstol": 1.0e-24,
                            "max_update": 2.0,
                            "handoff": {
                                "fallback": "none",
                                "require_gummel_convergence": False,
                                "gummel_max_iter": 0,
                                "newton_max_iter": 80,
                            },
                        },
                        "execute": False,
                    },
                    {
                        "name": "iv",
                        "kind": "iv",
                        "tdr": "pn2d_iv_des.tdr",
                        "cmd": "pn2d_iv_sdevice.cmd",
                        "plt": "pn2d_iv.plt",
                        "bias_column": "Anode OuterVoltage",
                        "current_column": "Anode TotalCurrent",
                        "vela_current_contact": "Cathode",
                        "execute": False,
                        "comparison": {
                            "candidate_column": "current_total_A_per_um",
                            "min_points": 2,
                        },
                    },
                ],
            }, indent=2) + "\n")
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
parser.add_argument("--compensated-doping-policy", default="reported")
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
    write_csv(out / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0, 0], [1, 1, 0], [2, 2, 0]])
    write_csv(out / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [[0, 0, 1, 2, "R.Si", "Si"]])
    write_csv(out / "contacts.csv", ["name", "node_ids", "region"], [["Anode", "0", "R.Si"], ["Cathode", "2", "R.Si"]])
    write_csv(out / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [[0, 0, 1e17], [1, 0, 1e17], [2, 1e17, 0]])
    write_csv(out / "fields" / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [[0, -0.4], [1, 0.0], [2, 0.4]])
    (out / "doping_metadata.json").write_text(json.dumps({"compensated_nodes": {"policy": args.compensated_doping_policy}}, indent=2) + "\\n")
    (out / "metadata.json").write_text(json.dumps({"vertex_count": 3, "region_count": 1, "dataset_count": 1}, indent=2) + "\\n")
    (out / "field_manifest.json").write_text(json.dumps({"fields": [{"name": "ElectrostaticPotential", "region": 0, "components": 1, "values": 3, "mapping_status": "complete"}]}, indent=2) + "\\n")
if args.inventory_json:
    Path(args.inventory_json).write_text(json.dumps({"vertex_count": 3, "region_count": 1, "dataset_count": 1}, indent=2) + "\\n")
""".strip()
                + "\n"
            )
            fake_runner = root / "fake_runner.py"
            fake_runner.write_text(
                """
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg = json.loads(Path(args.config).read_text())
out = Path(args.config).parent / cfg["output_csv"]
with out.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["bias_V", "current_total", "current_total_A_per_um"])
    writer.writerows([[0, 0, 0], [1.0, -1.25e-4, -1.25e-10]])
""".strip()
                + "\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "reference",
                    "--config",
                    str(config),
                    "--source-dir",
                    str(source),
                    "--output-dir",
                    str(output),
                    "--tdr-importer",
                    f"{sys.executable} {fake_importer}",
                    "--runner",
                    f"{sys.executable} {fake_runner}",
                ],
                check=True,
                cwd=REPO,
                env={**os.environ, "PYTHONIOENCODING": "utf-8"},
            )

            self.assertTrue((output / "reference_curves" / "pn2d_sentaurus2018_0v_reference.csv").is_file())
            self.assertTrue((output / "cmd" / "0v_summary.json").is_file())
            self.assertTrue((output / "sim_fields" / "0v" / "field_manifest.json").is_file())
            zero_deck = json.loads((output / "vela" / "simulation_0v.json").read_text())
            self.assertEqual(zero_deck["sweep"]["mode"], "equilibrium")
            self.assertEqual(zero_deck["sweep"]["start"], 0.0)
            self.assertEqual(zero_deck["sweep"]["stop"], 0.0)
            self.assertEqual(zero_deck["solver"]["method"], "gummel_newton")
            self.assertEqual(zero_deck["solver"]["max_iter"], 80)
            self.assertEqual(zero_deck["solver"]["reltol"], 1.0e-10)
            self.assertEqual(zero_deck["solver"]["abstol"], 1.0e-24)
            self.assertEqual(zero_deck["solver"]["max_update"], 2.0)
            self.assertEqual(zero_deck["solver"]["handoff"]["fallback"], "none")
            self.assertEqual(zero_deck["solver"]["handoff"]["gummel_max_iter"], 0)
            self.assertEqual(zero_deck["solver"]["handoff"]["newton_max_iter"], 80)
            manifest = json.loads((output / "reference_tcad_manifest.json").read_text())
            self.assertIn("0v execution disabled by reference config", manifest["warnings"])

    def test_compare_sentaurus_tdr_tdx_reports_mesh_and_field_parity(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_tdr_tdx_") as tmp:
            root = Path(tmp)
            export = root / "tdr_export"
            tdx = root / "tdx"
            reports = root / "reports"
            (export / "fields").mkdir(parents=True)
            tdx.mkdir()
            self._write_csv(export / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            self._write_csv(export / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
                [0, 0, 1, 2, "R.Si", "Si"],
            ])
            self._write_csv(export / "contacts.csv", ["name", "node_ids", "region"], [
                ["Anode", "0;2", "R.Si"],
            ])
            self._write_csv(export / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            self._write_csv(export / "fields" / "DopingConcentration_region0.csv", ["node_id", "component0"], [
                [0, -1.0e17],
                [1, 1.0e17],
                [2, -1.0e17],
            ])
            (tdx / "pn2d_msh.grd").write_text(
                """
DF-ISE text
Info {
  type = grid
  dimension = 2
  nb_vertices = 3
  nb_elements = 1
  nb_regions = 2
  regions = [ "R.Si" "Anode" ]
  materials = [ Silicon Contact ]
}
Data {
  Vertices (3) {
    0 0
    1 0
    0 1
  }
  Elements (2) {
    2 0 1 2
    1 0 2
  }
  Region ("R.Si") {
    material = Silicon
    Elements (1) { 0 }
  }
  Region ("Anode") {
    material = Contact
    Elements (1) { 1 }
  }
}
""".strip()
                + "\n"
            )
            (tdx / "pn2d_msh.dat").write_text(
                """
DF-ISE text
Info {
  type = dataset
  datasets = [ "DopingConcentration" ]
}
Data {
  Dataset ("DopingConcentration") {
    function = DopingConcentration
    type = scalar
    dimension = 1
    location = vertex
    validity = [ "R.Si" ]
    Values (3) {
      -1e17 1e17 -1e17
    }
  }
}
""".strip()
                + "\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_sentaurus_tdr_tdx.py"),
                    "--tdr-export",
                    str(export),
                    "--tdx-dir",
                    str(tdx),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "tdr_tdx_comparison.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["mesh"]["vertex_count"]["tdr"], 3)
            self.assertEqual(report["mesh"]["contact_element_counts"]["Anode"]["tdx"], 1)
            self.assertEqual(report["fields"]["DopingConcentration"]["max_abs_diff"], 0.0)

    def test_compare_sentaurus_tdr_tdx_reports_geometry_parity(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_geometry_") as tmp:
            root = Path(tmp)
            export = root / "tdr_export"
            tdx = root / "tdx"
            reports = root / "reports"
            self._write_geometry_fixture(export, tdx)

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_sentaurus_tdr_tdx.py"),
                    "--geometry-only",
                    "--tdr-export",
                    str(export),
                    "--tdx-dir",
                    str(tdx),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "tdr_tdx_comparison.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(sorted(report), ["boundaries", "elements", "nodes", "status"])
            self.assertTrue(report["nodes"]["count_match"])
            self.assertEqual(report["nodes"]["max_abs_diff"], 0.0)
            self.assertTrue(report["elements"]["count_match"])
            self.assertIsNone(report["elements"]["first_mismatch"])
            self.assertEqual(report["boundaries"]["contact_names"]["tdr"], ["Anode"])
            self.assertTrue(report["boundaries"]["contacts"]["Anode"]["node_set_match"])

    def test_compare_sentaurus_tdr_tdx_geometry_fails_on_integer_topology_mismatch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_geometry_int_fail_") as tmp:
            root = Path(tmp)
            export = root / "tdr_export"
            tdx = root / "tdx"
            reports = root / "reports"
            self._write_geometry_fixture(export, tdx, tdr_triangle=[0, 1, 3])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_sentaurus_tdr_tdx.py"),
                    "--geometry-only",
                    "--tdr-export",
                    str(export),
                    "--tdx-dir",
                    str(tdx),
                    "--output-dir",
                    str(reports),
                ],
                cwd=REPO,
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((reports / "tdr_tdx_comparison.json").read_text())
            self.assertEqual(report["status"], "fail")
            self.assertEqual(report["elements"]["first_mismatch"]["tdr_nodes"], [0, 1, 3])
            self.assertEqual(report["elements"]["first_mismatch"]["grd_nodes"], [0, 1, 2])

    def test_compare_sentaurus_tdr_tdx_geometry_fails_on_coordinate_error_above_epsilon(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_geometry_float_fail_") as tmp:
            root = Path(tmp)
            export = root / "tdr_export"
            tdx = root / "tdx"
            reports = root / "reports"
            self._write_geometry_fixture(export, tdx, tdr_node1_x=1.0 + 1.0e-12)

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_sentaurus_tdr_tdx.py"),
                    "--geometry-only",
                    "--tdr-export",
                    str(export),
                    "--tdx-dir",
                    str(tdx),
                    "--output-dir",
                    str(reports),
                ],
                cwd=REPO,
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((reports / "tdr_tdx_comparison.json").read_text())
            self.assertEqual(report["status"], "fail")
            self.assertEqual(report["nodes"]["worst"]["node_id"], 1)
            self.assertEqual(report["nodes"]["worst"]["component"], "x_um")

    def test_compare_sentaurus_fields_reports_selected_field_parity(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_field_compare_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            candidate = root / "candidate"
            reports = root / "reports"
            reference.mkdir()
            candidate.mkdir()
            self._write_csv(reference / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [0, -0.4],
                [1, 0.0],
                [2, 0.4],
            ])
            self._write_csv(candidate / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [0, -0.4000000001],
                [1, 0.0],
                [2, 0.4000000001],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_sentaurus_fields.py"),
                    "--reference-fields",
                    str(reference),
                    "--candidate-fields",
                    str(candidate),
                    "--fields",
                    "ElectrostaticPotential",
                    "--output-json",
                    str(reports / "field_comparison.json"),
                    "--max-abs-diff",
                    "1e-8",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "field_comparison.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["fields"]["ElectrostaticPotential"]["points_compared"], 3)
            self.assertLess(report["fields"]["ElectrostaticPotential"]["max_abs_diff"], 1.0e-8)

    def test_compare_pn2d_0v_state_reports_field_formula_and_terminal_parity(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_state_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            vela = reference / "vela"
            reports = root / "reports"
            fields.mkdir(parents=True)
            vela.mkdir(parents=True)
            self._write_csv(reference / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 2.0, 0.0],
            ])
            self._write_csv(reference / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e10],
                [1, 1.0e10, 0.0],
                [2, 1.0e10, 0.0],
            ])
            self._write_csv(reference / "contacts.csv", ["name", "node_ids", "region"], [
                ["Anode", "0", "R.Si"],
                ["Cathode", "2", "R.Si"],
            ])
            self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [0, -0.025851999786435],
                [1, 0.0],
                [2, 0.025851999786435],
            ])
            self._write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
                [2, 0.0],
            ])
            self._write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
                [2, 0.0],
            ])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [
                [0, 3678794411.714423],
                [1, 10000000000.0],
                [2, 27182818284.59045],
            ])
            self._write_csv(fields / "hDensity_region0.csv", ["node_id", "component0"], [
                [0, 27182818284.59045],
                [1, 10000000000.0],
                [2, 3678794411.714423],
            ])
            self._write_csv(fields / "srhRecombination_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
                [2, 0.0],
            ])
            (vela / "simulation_0v.json").write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "output_csv": "pn2d_sentaurus2018_0v.csv",
                "scaling": {"mode": "unit_scaling"},
                "node_doping_file": "doping.csv",
                "contacts": [{"name": "Anode", "bias": 0.0}, {"name": "Cathode", "bias": 0.0}],
                "solver": {
                    "method": "gummel_newton",
                    "recombination": ["srh"],
                    "bandgap_narrowing": "none",
                    "handoff": {"fallback": "none", "gummel_max_iter": 0, "newton_max_iter": 5},
                },
                "sweep": {
                    "mode": "equilibrium",
                    "contact": "Anode",
                    "current_contact": "Anode",
                    "start": 0.0,
                    "stop": 0.0,
                    "step": 0.0,
                    "write_vtk": False,
                },
            }, indent=2) + "\n")
            vtk = root / "vela_0v.vtk"
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
2 0 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
-0.025851999786435 0.0 0.025851999786435
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0.0 0.0 0.0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0.0 0.0 0.0
SCALARS Electrons float 1
LOOKUP_TABLE default
3.678794411714423e15 1.0e16 2.718281828459045e16
SCALARS Holes float 1
LOOKUP_TABLE default
2.718281828459045e16 1.0e16 3.678794411714423e15
""".lstrip()
            )
            terminal = root / "terminal.csv"
            self._write_csv(terminal, [
                "mode", "bias_contact", "bias_V", "current_contact",
                "current_total", "converged", "solver_method",
                "newton_iterations", "handoff_stage", "current_total_A_per_um",
            ], [
                ["iv", "Anode", 0.0, "Anode", 1.0e-18, 1, "gummel_newton", 3, "newton", 1.0e-24],
                ["iv", "Anode", 0.0, "Cathode", -1.0e-18, 1, "gummel_newton", 3, "newton", -1.0e-24],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_0v_state.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-vtk",
                    str(vtk),
                    "--existing-terminal-csv",
                    str(terminal),
                    "--output-dir",
                    str(reports),
                    "--ni-cm3",
                    "1e10",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_state_comparison.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["field_stats"]["ElectrostaticPotential"]["points_compared"], 3)
            self.assertLess(report["field_stats"]["ElectrostaticPotential"]["max_abs_diff"], 1.0e-12)
            self.assertLess(report["field_stats"]["eDensity_formula"]["max_rel_diff"], 1.0e-12)
            self.assertLess(report["field_stats"]["srhRecombination"]["max_abs_diff"], 1.0e-6)
            self.assertEqual(report["diagnostic_matrix"]["ni_cm3"], [1.0e10])
            self.assertEqual(report["diagnostic_matrix"]["bandgap_narrowing"], ["none", "slotboom"])
            self.assertEqual(report["terminal_currents"]["status"], "pass")
            self.assertLessEqual(abs(report["terminal_currents"]["sum_A_per_um"]), 1.0e-21)

    def test_compare_pn2d_0v_state_fails_when_terminal_currents_are_not_balanced(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_state_terminal_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            reports = root / "reports"
            fields.mkdir(parents=True)
            self._write_csv(reference / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0]])
            self._write_csv(reference / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [[0, 1.0e10, 0.0]])
            for name in [
                "ElectrostaticPotential",
                "eQuasiFermiPotential",
                "hQuasiFermiPotential",
                "srhRecombination",
            ]:
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [[0, 0.0]])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [[0, 1.0e10]])
            self._write_csv(fields / "hDensity_region0.csv", ["node_id", "component0"], [[0, 1.0e10]])
            vtk = root / "minimal.vtk"
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 1 float
0 0 0
POINT_DATA 1
SCALARS Potential float 1
LOOKUP_TABLE default
0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0
""".lstrip()
            )
            terminal = root / "terminal.csv"
            self._write_csv(terminal, ["current_contact", "current_total_A_per_um"], [
                ["Anode", 4.0e-21],
                ["Cathode", 2.0e-27],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_0v_state.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-vtk",
                    str(vtk),
                    "--existing-terminal-csv",
                    str(terminal),
                    "--output-dir",
                    str(reports),
                ],
                cwd=REPO,
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((reports / "pn2d_0v_state_comparison.json").read_text())
            self.assertEqual(report["status"], "fail")
            self.assertEqual(report["terminal_currents"]["status"], "fail")
            self.assertGreater(report["terminal_currents"]["pair_balance_relative"], 0.95)

    def test_compare_pn2d_0v_state_fails_with_missing_required_field(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_state_missing_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            vela = reference / "vela"
            reports = root / "reports"
            fields.mkdir(parents=True)
            vela.mkdir(parents=True)
            self._write_csv(reference / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0]])
            self._write_csv(reference / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [[0, 0.0, 1.0e10]])
            self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [[0, 0.0]])
            self._write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [[0, 0.0]])
            (vela / "simulation_0v.json").write_text("{}\n")
            vtk = root / "minimal.vtk"
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 1 float
0 0 0
POINT_DATA 1
SCALARS Potential float 1
LOOKUP_TABLE default
0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0
""".lstrip()
            )
            terminal = root / "terminal.csv"
            self._write_csv(terminal, ["current_contact", "current_total_A_per_um"], [["Anode", 0.0]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_0v_state.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-vtk",
                    str(vtk),
                    "--existing-terminal-csv",
                    str(terminal),
                    "--output-dir",
                    str(reports),
                ],
                cwd=REPO,
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((reports / "pn2d_0v_state_comparison.json").read_text())
            self.assertEqual(report["status"], "fail")
            self.assertIn("sentaurus:hQuasiFermiPotential", report["missing_fields"])
            self.assertIn("vtk:HoleQuasiFermi", report["missing_fields"])

    def test_diagnose_pn2d_0v_field_conventions_identifies_potential_sign_and_density_formula(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_field_conventions_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            vela = reference / "vela"
            reports = root / "reports"
            fields.mkdir(parents=True)
            vela.mkdir(parents=True)
            vt = 8.617333262145e-5 * 300.0
            ni = 1.0e10
            sentaurus_psi = [0.3, 0.1, -0.1]
            vela_psi = [-0.2, 0.0, 0.2]
            self._write_csv(reference / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e10],
                [1, 0.0, 0.0],
                [2, 1.0e10, 0.0],
            ])
            self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [idx, value] for idx, value in enumerate(sentaurus_psi)
            ])
            self._write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0], [1, 0.0], [2, 0.0]
            ])
            self._write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0], [1, 0.0], [2, 0.0]
            ])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [
                [idx, ni * math.exp(value / vt)] for idx, value in enumerate(sentaurus_psi)
            ])
            self._write_csv(fields / "hDensity_region0.csv", ["node_id", "component0"], [
                [idx, ni * math.exp(-value / vt)] for idx, value in enumerate(sentaurus_psi)
            ])
            (vela / "simulation_0v.json").write_text(json.dumps({
                "solver": {"bandgap_narrowing": "none"}
            }, indent=2) + "\n")
            vtk = root / "vela.vtk"
            vtk.write_text(
                f"""
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
2 0 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
{vela_psi[0]} {vela_psi[1]} {vela_psi[2]}
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_field_conventions.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--ni-cm3",
                    "1e10",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_field_conventions.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["potential_convention"]["best"], "opposite_sign_with_offset")
            self.assertLess(report["potential_convention"]["candidates"]["opposite_sign_with_offset"]["max_abs_diff"], 1.0e-12)
            self.assertEqual(report["carrier_formula"]["electron"]["best"]["formula"], "psi_minus_qf")
            self.assertEqual(report["carrier_formula"]["hole"]["best"]["formula"], "qf_minus_psi")
            self.assertLess(report["carrier_formula"]["electron"]["best"]["stats"]["max_rel_diff"], 1.0e-12)
            self.assertLess(report["carrier_formula"]["hole"]["best"]["stats"]["max_rel_diff"], 1.0e-12)

    def test_diagnose_pn2d_0v_field_mapping_detects_coordinate_ordering(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_field_mapping_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            reports = root / "reports"
            fields.mkdir(parents=True)
            self._write_csv(reference / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 2.0, 0.0],
            ])
            self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [0, 10.0],
                [1, 20.0],
                [2, 30.0],
            ])
            vtk = root / "permuted.vtk"
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
2 0 0
0 0 0
1 0 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
30 10 20
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_field_mapping.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--fields",
                    "ElectrostaticPotential:Potential",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_field_mapping.json").read_text())
            field = report["fields"]["ElectrostaticPotential"]
            self.assertEqual(field["best_pairing"], "nearest_coordinate")
            self.assertEqual(field["direct_node_id"]["max_abs_diff"], 20.0)
            self.assertLess(field["nearest_coordinate"]["max_abs_diff"], 1.0e-12)
            self.assertEqual(field["nearest_coordinate"]["matched_points"], 3)
            self.assertEqual(report["coordinate_alignment"]["best_scale"], "um")

    def test_diagnose_pn2d_0v_field_mapping_identifies_density_units(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_density_units_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            reports = root / "reports"
            fields.mkdir(parents=True)
            self._write_csv(reference / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 2.0, 0.0],
            ])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [
                [0, 1.0e10],
                [1, 2.0e10],
                [2, 4.0e10],
            ])
            vtk = root / "density_si.vtk"
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
2 0 0
POINT_DATA 3
SCALARS Electrons float 1
LOOKUP_TABLE default
1.0e16 2.0e16 4.0e16
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_field_mapping.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--fields",
                    "eDensity:Electrons",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_field_mapping.json").read_text())
            field = report["fields"]["eDensity"]
            self.assertEqual(field["best_unit_convention"], "vtk_m3_to_cm3")
            self.assertEqual(field["unit_classification"], "vtk_m3_reference_cm3")
            self.assertEqual(field["unit_conventions"]["raw"]["median_candidate_over_reference"], 1.0e6)
            self.assertLess(field["unit_conventions"]["vtk_m3_to_cm3"]["max_abs_diff"], 1.0e-12)

    def test_compare_pn2d_0v_electric_field_derives_v_per_cm_from_potential(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_electric_field_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            reports = root / "reports"
            fields.mkdir(parents=True)
            (reference / "vela").mkdir()
            self._write_csv(reference / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0e6, 0.0],
                [2, 0.0, 1.0e6],
            ])
            (reference / "vela" / "mesh.json").write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                ],
                "contacts": [
                    {"name": "Anode", "node_ids": [0]},
                    {"name": "Cathode", "node_ids": [1]},
                ],
            }))
            self._write_csv(fields / "ElectricField_region0.csv", ["node_id", "component0"], [
                [0, 100.0],
                [1, 100.0],
                [2, 100.0],
            ])
            vtk = root / "vela_state.vtk"
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
0 1 0
POLYGONS 1 4
3 0 1 2
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 10000 0
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_0v_electric_field.py"),
                    "--reference-root",
                    str(reference),
                    "--vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_electric_field_comparison.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["source"], "derived_from_potential_gradient")
            self.assertEqual(report["unit_selection"]["best_unit"], "V_per_cm")
            self.assertLess(report["stats"]["best_unit"]["all"]["max_abs_diff"], 1.0e-12)
            self.assertEqual(report["stats"]["best_unit"]["contact:Anode"]["points"], 1)
            with (reports / "pn2d_0v_electric_field_nodes.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["unit"], "V_per_cm")

    def test_compare_pn2d_0v_electric_field_symmetry_reports_solver_asymmetry(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_electric_field_symmetry_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            reports = root / "reports"
            field_report = reference / "reports" / "0v_electric_field"
            field_report.mkdir(parents=True)
            (reference / "vela").mkdir()
            self._write_csv(reference / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 2.0, 0.0],
                [2, 0.0, 1.0],
                [3, 2.0, 1.0],
            ])
            (reference / "vela" / "mesh.json").write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 2.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 2.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                    {"id": 1, "node_ids": [1, 3, 2]},
                ],
            }))
            self._write_csv(
                field_report / "pn2d_0v_electric_field_nodes.csv",
                ["node_id", "unit", "sentaurus_E", "vela_E", "abs_diff", "rel_diff"],
                [
                    [0, "V_per_cm", 10.0, 10.0, 0.0, 0.0],
                    [1, "V_per_cm", 10.0, 30.0, 20.0, 2.0],
                    [2, "V_per_cm", 20.0, 20.0, 0.0, 0.0],
                    [3, "V_per_cm", 20.0, 80.0, 60.0, 3.0],
                ],
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_0v_electric_field_symmetry.py"),
                    "--reference-root",
                    str(reference),
                    "--output-dir",
                    str(reports),
                    "--axis-x-um",
                    "1.0",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_electric_field_symmetry.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["axis_x_um"], 1.0)
            self.assertEqual(report["pairing"]["pairs"], 2)
            self.assertEqual(report["solvers"]["sentaurus"]["max_abs_diff_V_per_cm"], 0.0)
            self.assertEqual(report["solvers"]["sentaurus"]["classification"], "symmetric")
            self.assertEqual(report["solvers"]["vela"]["max_abs_diff_V_per_cm"], 60.0)
            self.assertEqual(report["solvers"]["vela"]["classification"], "asymmetric")
            self.assertTrue((reports / "pn2d_0v_sentaurus_electric_field_symmetry.png").is_file())
            self.assertTrue((reports / "pn2d_0v_vela_electric_field_symmetry.png").is_file())

    def test_diagnose_pn2d_0v_density_decomposition_classifies_state_difference(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_density_decomposition_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            reports = root / "reports"
            fields.mkdir(parents=True)
            ni = 1.0e10
            vt = 8.617333262145e-5 * 300.0
            sentaurus_psi = [0.10, -0.10]
            vela_psi = [0.0, 0.0]
            self._write_csv(reference / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 0.0, 1.0e17],
            ])
            self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [0, sentaurus_psi[0]],
                [1, sentaurus_psi[1]],
            ])
            self._write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
            ])
            self._write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
            ])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [
                [idx, ni * math.exp(value / vt)] for idx, value in enumerate(sentaurus_psi)
            ])
            self._write_csv(fields / "hDensity_region0.csv", ["node_id", "component0"], [
                [idx, ni * math.exp(-value / vt)] for idx, value in enumerate(sentaurus_psi)
            ])
            vtk = root / "vela_state.vtk"
            vtk.write_text(
                f"""
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0 0 0
1 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
{vela_psi[0]} {vela_psi[1]}
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
{ni * math.exp(vela_psi[0] / vt) * 1.0e6} {ni * math.exp(vela_psi[1] / vt) * 1.0e6}
SCALARS Holes float 1
LOOKUP_TABLE default
{ni * math.exp(-vela_psi[0] / vt) * 1.0e6} {ni * math.exp(-vela_psi[1] / vt) * 1.0e6}
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_density_decomposition.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--ni-cm3",
                    str(ni),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_density_decomposition.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["classification"], "vela_state_differs_from_sentaurus")
            self.assertLess(report["vela_self_consistency"]["electron"]["max_rel_diff"], 1.0e-12)
            self.assertLess(report["vela_self_consistency"]["hole"]["max_rel_diff"], 1.0e-12)
            self.assertLess(report["sentaurus_self_consistency"]["electron"]["max_rel_diff"], 1.0e-12)
            self.assertGreater(report["sentaurus_vs_vela_density"]["electron_vtk_cm3"]["max_rel_diff"], 1.0)

    def test_diagnose_pn2d_0v_ni_bgn_probe_ranks_material_ni_candidate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_0v_ni_probe_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            fields = reference / "sim_fields" / "0v" / "fields"
            vela = reference / "vela"
            reports = root / "reports"
            fields.mkdir(parents=True)
            vela.mkdir(parents=True)
            target_ni = 1.6556207295e10
            self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
            ])
            self._write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
            ])
            self._write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 0.0],
            ])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [
                [0, target_ni],
                [1, target_ni],
            ])
            self._write_csv(fields / "hDensity_region0.csv", ["node_id", "component0"], [
                [0, target_ni],
                [1, target_ni],
            ])
            self._write_csv(reference / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 0.0, 1.0e17],
            ])
            (vela / "simulation_0v.json").write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "solver": {"bandgap_narrowing": "none"},
                "sweep": {"mode": "iv", "start": 0.0, "stop": 0.0, "step": 1.0},
            }, indent=2) + "\n")
            fake_runner = root / "fake_runner.py"
            fake_runner.write_text(
                """
import json
import sys
from pathlib import Path

config = Path(sys.argv[sys.argv.index("--config") + 1])
deck = json.loads(config.read_text())
materials = json.loads(Path(deck["materials_file"]).read_text())
ni = float(materials["materials"][0]["ni"])
prefix = Path(deck["sweep"]["vtk_prefix"])
vtk = prefix.with_name(prefix.name + "_0.vtk")
vtk.parent.mkdir(parents=True, exist_ok=True)
vtk.write_text(f'''# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0 0 0
1 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
{ni * 1e6} {ni * 1e6}
SCALARS Holes float 1
LOOKUP_TABLE default
{ni * 1e6} {ni * 1e6}
''')
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_ni_bgn_probe.py"),
                    "--reference-root",
                    str(reference),
                    "--runner",
                    f"{sys.executable} {fake_runner}",
                    "--output-dir",
                    str(reports),
                    "--ni-cm3",
                    "1e10",
                    "--ni-cm3",
                    str(target_ni),
                    "--bgn",
                    "none",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_ni_bgn_probe.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertAlmostEqual(report["best_candidate"]["ni_cm3"], target_ni)
            self.assertEqual(report["best_candidate"]["bandgap_narrowing"], "none")
            self.assertLess(report["best_candidate"]["metrics"]["density"]["max_log10_error"], 1.0e-12)

    def test_diagnose_pn2d_0v_current_balance_classifies_sign_convention(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_sign_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 10.0, 3.0, 7.0, 7.0, 13.0],
                ["Cathode", -4.0, -9.0, 5.0, 5.0, -13.0],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 7.0], ["Cathode", 1, 5.0]])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["classification"], "total_current_sign_convention")
            self.assertTrue(report["root_cause_flags"]["total_current_sign_convention"])
            self.assertEqual(report["terminal_balance"]["conventions"]["electron_plus_hole"]["status"], "pass")
            self.assertEqual(report["terminal_balance"]["conventions"]["electron_minus_hole"]["status"], "fail")
            self.assertEqual(report["mesh_reference"]["sentaurus_contact_boundary_elements"]["Anode"], 1)
            self.assertEqual(report["mesh_reference"]["sentaurus_contact_boundary_elements"]["Cathode"], 1)
            self.assertEqual(report["contact_edges"]["flux_link_count_by_contact"]["Anode"], 1)

    def test_diagnose_pn2d_0v_current_balance_classifies_edge_aggregation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_aggregation_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 10.0, 3.0, 7.0, 7.0, 13.0],
                ["Cathode", -10.0, -3.0, -7.0, -7.0, -13.0],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 6.0], ["Cathode", 1, -7.0]])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["classification"], "contact_current_aggregation")
            self.assertTrue(report["root_cause_flags"]["contact_current_aggregation"])
            self.assertEqual(report["contact_edges"]["aggregation"]["Anode"]["status"], "fail")

    def test_diagnose_pn2d_0v_current_balance_classifies_edge_coverage(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_coverage_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 10.0, 3.0, 7.0, 7.0, 13.0],
                ["Cathode", -10.0, -3.0, -7.0, -7.0, -13.0],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 7.0]])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["classification"], "contact_edge_coverage")
            self.assertTrue(report["root_cause_flags"]["contact_edge_coverage"])
            self.assertIn("Cathode", report["contact_edges"]["missing_contacts"])

    def test_diagnose_pn2d_0v_current_balance_require_balanced_passes_on_absolute_floor(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_require_pass_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 2.0e-30, 1.0e-30, 1.0e-30, 1.0e-30, 3.0e-30],
                ["Cathode", -1.0e-30, -1.0e-30, 0.0, 0.0, -2.0e-30],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 1.0e-30], ["Cathode", 1, 0.0]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--require-balanced",
                ],
                cwd=REPO,
            )

            self.assertEqual(result.returncode, 0)
            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["classification"], "balanced")
            summary = report["conservation_summary"]
            self.assertLessEqual(summary["abs_pair_sum_A_per_um"], 1.0e-24)
            self.assertTrue(summary["absolute_floor_gate_pass"])

    def test_diagnose_pn2d_0v_current_balance_require_balanced_fails_for_sign_convention(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_require_sign_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 10.0, 3.0, 7.0, 7.0, 13.0],
                ["Cathode", -4.0, -9.0, 5.0, 5.0, -13.0],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 7.0], ["Cathode", 1, 5.0]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--require-balanced",
                ],
                cwd=REPO,
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["classification"], "total_current_sign_convention")

    def test_diagnose_pn2d_0v_current_balance_require_balanced_fails_for_aggregation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_require_aggregation_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 10.0, 3.0, 7.0, 7.0, 13.0],
                ["Cathode", -10.0, -3.0, -7.0, -7.0, -13.0],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 6.0], ["Cathode", 1, -7.0]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--require-balanced",
                ],
                cwd=REPO,
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["classification"], "contact_current_aggregation")

    def test_diagnose_pn2d_0v_current_balance_require_balanced_fails_for_edge_coverage(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_require_coverage_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 10.0, 3.0, 7.0, 7.0, 13.0],
                ["Cathode", -10.0, -3.0, -7.0, -7.0, -13.0],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 7.0]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                    "--require-balanced",
                ],
                cwd=REPO,
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["classification"], "contact_edge_coverage")

    def test_diagnose_pn2d_0v_current_balance_reports_sentaurus_current_parity(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_sentaurus_parity_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            self._write_sentaurus_0v_reference_files(reference)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 1.0e-33, 0.0, 1.0e-33, 1.0e-33, 1.0e-33],
                ["Cathode", 0.0, 8.0e-33, -8.0e-33, -8.0e-33, 8.0e-33],
            ])
            self._write_contact_edge_csv(edges, [["Anode", 0, 1.0e-33], ["Cathode", 1, -8.0e-33]])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            sentaurus = report["sentaurus_current_reference"]
            self.assertEqual(sentaurus["final_coupled"]["source"], "plt")
            self.assertAlmostEqual(sentaurus["final_coupled"]["contacts"]["Anode"]["total_current"], -7.17389811693691e-25)
            self.assertAlmostEqual(sentaurus["final_coupled"]["contacts"]["Cathode"]["total_current"], 7.17389811693687e-25)
            parity = report["sentaurus_current_parity"]
            self.assertEqual(parity["status"], "mismatch")
            self.assertEqual(parity["unit_hypothesis"], "vela_current_total_A_per_um_vs_sentaurus_current")
            self.assertGreater(parity["by_contact"]["Anode"]["abs_ratio_sentaurus_to_vela"], 1.0e8)
            self.assertIn("Anode", parity["minority_component_zero_contacts"])
            self.assertIn("Cathode", parity["minority_component_zero_contacts"])
            self.assertEqual(parity["classification"], "sentaurus_abs_current_mismatch")
            self.assertEqual(
                report["contact_edges"]["component_sum_A_per_um_by_contact"]["Anode"]["electron_minus_hole"],
                1.0e-33,
            )
            self.assertEqual(
                report["contact_edges"]["component_sum_A_per_um_by_contact"]["Cathode"]["electron_minus_hole"],
                -8.0e-33,
            )

    def test_diagnose_pn2d_0v_current_balance_explains_zero_qf_edge_component(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_balance_zero_qf_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_current_balance_common_fixture(root)
            terminal = root / "terminal_balance.csv"
            edges = root / "contact_edges.csv"
            reports = root / "reports"
            self._write_terminal_balance_csv(terminal, [
                ["Anode", 1.0e-33, 0.0, 1.0e-33, 1.0e-33, 1.0e-33],
                ["Cathode", 0.0, 1.0e-33, -1.0e-33, -1.0e-33, 1.0e-33],
            ])
            self._write_csv(edges, [
                "current_contact",
                "edge_id",
                "node0",
                "node1",
                "current_total_A_per_um",
                "psi0",
                "psi1",
                "phin0",
                "phin1",
                "phip0",
                "phip1",
                "n0",
                "n1",
                "p0",
                "p1",
                "ni0",
                "ni1",
                "mun",
                "mup",
                "electron_branch",
                "hole_branch",
                "electron_continuity_flux",
                "hole_continuity_flux",
                "current_electron",
                "current_hole",
            ], [
                ["Anode", 0, 0, 1, 1.0e-33, 0, 0, 1.0e-16, 0.0, 1.0e-19, 0.0,
                 1.0e9, 1.0e9, 1.0e23, 1.0e23, 1.0e16, 1.0e16, 0.135, 0.048,
                 "quasi_fermi", "quasi_fermi", 1.0e-2, 0.0, 1.0e-27, 0.0],
                ["Cathode", 1, 0, 1, -1.0e-33, 0, 0, 1.0e-19, 0.0, 1.0e-16, 0.0,
                 1.0e23, 1.0e23, 1.0e9, 1.0e9, 1.0e16, 1.0e16, 0.135, 0.048,
                 "quasi_fermi", "quasi_fermi", 0.0, 1.0e-2, 0.0, 1.0e-27],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(reference),
                    "--existing-terminal-balance-csv",
                    str(terminal),
                    "--existing-contact-edge-csv",
                    str(edges),
                    "--existing-vtk",
                    str(vtk),
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            zero_components = report["contact_edges"]["zero_component_diagnostics"]
            self.assertEqual(
                zero_components["by_contact"]["Anode"]["hole"]["dominant_reason"],
                "qf_delta_below_exp_resolution",
            )
            self.assertEqual(
                zero_components["by_contact"]["Cathode"]["electron"]["dominant_reason"],
                "qf_delta_below_exp_resolution",
            )
            self.assertGreater(zero_components["by_contact"]["Anode"]["hole"]["count"], 0)

    def test_probe_pn2d_0v_newton_residual_current_ranks_relaxed_candidate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_newton_residual_probe_") as tmp:
            root = Path(tmp)
            reference, _vtk = self._write_current_balance_common_fixture(root)
            self._write_sentaurus_0v_reference_files(reference)
            fake_runner = root / "fake_runner.py"
            self._write_newton_residual_fake_runner(fake_runner)
            reports = root / "reports"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "probe_pn2d_0v_newton_residual_current.py"),
                    "--reference-root",
                    str(reference),
                    "--runner",
                    f"{sys.executable} {fake_runner}",
                    "--output-dir",
                    str(reports),
                    "--reltol-candidate",
                    "1e-10",
                    "--reltol-candidate",
                    "1e-2",
                    "--abstol-candidate",
                    "1e-24",
                    "--max-iter-candidate",
                    "80",
                    "--max-iter-candidate",
                    "1",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_newton_residual_probe.json").read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(len(report["candidates"]), 4)
            self.assertEqual(report["best_candidate"]["solver"]["reltol"], 1.0e-2)
            self.assertEqual(report["best_candidate"]["solver"]["max_iter"], 1)
            self.assertLess(report["best_candidate"]["sentaurus_total_abs_ratio_median_log10_error"], 0.1)
            strict = next(item for item in report["candidates"] if item["solver"]["reltol"] == 1.0e-10)
            self.assertGreater(strict["sentaurus_total_abs_ratio_median_log10_error"], 7.0)

    def test_probe_pn2d_0v_newton_residual_current_records_failed_candidates(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_newton_residual_failed_") as tmp:
            root = Path(tmp)
            reference, _vtk = self._write_current_balance_common_fixture(root)
            self._write_sentaurus_0v_reference_files(reference)
            fake_runner = root / "fake_runner.py"
            self._write_newton_residual_fake_runner(fake_runner)
            with fake_runner.open("a") as handle:
                handle.write(
                    "\nif int(solver.get('max_iter', 80)) == 3:\n"
                    "    raise SystemExit(7)\n"
                )
            reports = root / "reports"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "probe_pn2d_0v_newton_residual_current.py"),
                    "--reference-root",
                    str(reference),
                    "--runner",
                    f"{sys.executable} {fake_runner}",
                    "--output-dir",
                    str(reports),
                    "--reltol-candidate",
                    "1e-10",
                    "--abstol-candidate",
                    "1e-24",
                    "--max-iter-candidate",
                    "80",
                    "--max-iter-candidate",
                    "3",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_newton_residual_probe.json").read_text())
            failed = next(item for item in report["candidates"] if item["solver"]["max_iter"] == 3)
            self.assertEqual(failed["status"], "runner_failed")
            self.assertEqual(failed["sentaurus_total_abs_ratio_median_log10_error"], 300.0)
            self.assertEqual(report["best_candidate"]["solver"]["max_iter"], 80)

    def test_probe_pn2d_0v_newton_residual_current_prefers_balanced_candidate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_newton_residual_balance_rank_") as tmp:
            root = Path(tmp)
            reference, _vtk = self._write_current_balance_common_fixture(root)
            self._write_sentaurus_0v_reference_files(reference)
            fake_runner = root / "fake_runner.py"
            self._write_newton_residual_unbalanced_fake_runner(fake_runner)
            reports = root / "reports"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "probe_pn2d_0v_newton_residual_current.py"),
                    "--reference-root",
                    str(reference),
                    "--runner",
                    f"{sys.executable} {fake_runner}",
                    "--output-dir",
                    str(reports),
                    "--reltol-candidate",
                    "1e-10",
                    "--reltol-candidate",
                    "1e-2",
                    "--abstol-candidate",
                    "1e-24",
                    "--max-iter-candidate",
                    "80",
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_newton_residual_probe.json").read_text())
            self.assertEqual(report["best_candidate"]["terminal_balance"]["status"], "pass")
            self.assertEqual(report["best_candidate"]["solver"]["reltol"], 1.0e-10)

    def _write_current_balance_common_fixture(self, root: Path) -> tuple[Path, Path]:
        reference = root / "reference"
        (reference / "vela").mkdir(parents=True)
        (reference / "vela" / "simulation_0v.json").write_text("{}\n")
        (reference / "tdr_inventory").mkdir(parents=True)
        (reference / "tdr_inventory" / "mesh.json").write_text(json.dumps({
            "geometry": {
                "regions": [
                    {"name": "R.Si", "type": 0, "triangles": [[0, 1, 2], [1, 2, 3]]},
                    {"name": "Anode", "type": 1, "edges": [[0, 2]]},
                    {"name": "Cathode", "type": 1, "edges": [[1, 3]]},
                ]
            }
        }, indent=2) + "\n")
        (reference / "vela" / "mesh.json").write_text(json.dumps({
            "triangles": [
                {"id": 0, "node_ids": [0, 1, 2]},
                {"id": 1, "node_ids": [1, 2, 3]},
            ],
            "contacts": [
                {"name": "Anode", "node_ids": [0, 2]},
                {"name": "Cathode", "node_ids": [1, 3]},
            ],
        }, indent=2) + "\n")
        vtk = root / "balance.vtk"
        vtk.write_text(
            """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0 0 0
1 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0
""".lstrip()
        )
        return reference, vtk

    def _write_newton_residual_fake_runner(self, path: Path) -> None:
        path.write_text(
            """
import csv
import json
import sys
from pathlib import Path

cfg = json.loads(Path(sys.argv[sys.argv.index("--config") + 1]).read_text())
solver = cfg.get("solver", {})
reltol = float(solver.get("reltol", 1e-10))
max_iter = int(solver.get("max_iter", 80))
diag = cfg.get("sweep", {}).get("diagnostics", {})
terminal = Path(diag["terminal_balance"]["csv_file"])
edges = Path(diag["contact_edge"]["csv_file"])
terminal.parent.mkdir(parents=True, exist_ok=True)
edges.parent.mkdir(parents=True, exist_ok=True)

if reltol >= 1e-2 and max_iter <= 1:
    current = 7.1738981169369e-25
else:
    current = 1.0e-33

with terminal.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "point_index", "bias_V", "bias_contact", "contact",
        "current_electron", "current_hole", "electron_minus_hole",
        "current_total", "electron_plus_hole",
        "current_electron_A_per_um", "current_hole_A_per_um",
        "electron_minus_hole_A_per_um", "current_total_A_per_um",
        "electron_plus_hole_A_per_um", "converged", "solver_method",
        "gummel_iterations", "newton_iterations", "handoff_stage",
    ])
    writer.writerow([0, 0, "Anode", "Anode", 0, -current, current, current, -current,
                     0, -current, current, current, -current, 1, "gummel_newton", 0, max_iter, "newton"])
    writer.writerow([0, 0, "Anode", "Cathode", 0, current, -current, -current, current,
                     0, current, -current, -current, current, 1, "gummel_newton", 0, max_iter, "newton"])

with edges.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "current_contact", "edge_id", "node0", "node1", "current_total_A_per_um",
        "psi0", "psi1", "phin0", "phin1", "phip0", "phip1",
        "n0", "n1", "p0", "p1", "ni0", "ni1", "mun", "mup",
        "electron_branch", "hole_branch", "electron_continuity_flux",
        "hole_continuity_flux", "current_electron", "current_hole",
    ])
    writer.writerow(["Anode", 0, 0, 1, current, 0, 0, 0, 0, 0, 0,
                     1, 1, 1, 1, 1, 1, 0.135, 0.048, "quasi_fermi", "quasi_fermi", 0, -current, 0, -current])
    writer.writerow(["Cathode", 1, 0, 1, -current, 0, 0, 0, 0, 0, 0,
                     1, 1, 1, 1, 1, 1, 0.135, 0.048, "quasi_fermi", "quasi_fermi", 0, current, 0, current])

vtk_prefix = Path(cfg.get("sweep", {}).get("vtk_prefix", terminal.parent / "probe"))
vtk = vtk_prefix.with_name(vtk_prefix.name + "_0000.vtk")
vtk.write_text(\"\"\"# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0 0 0
1 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0
\"\"\")
""".lstrip()
        )

    def _write_newton_residual_unbalanced_fake_runner(self, path: Path) -> None:
        path.write_text(
            """
import csv
import json
import sys
from pathlib import Path

cfg = json.loads(Path(sys.argv[sys.argv.index("--config") + 1]).read_text())
solver = cfg.get("solver", {})
reltol = float(solver.get("reltol", 1e-10))
diag = cfg.get("sweep", {}).get("diagnostics", {})
terminal = Path(diag["terminal_balance"]["csv_file"])
edges = Path(diag["contact_edge"]["csv_file"])
terminal.parent.mkdir(parents=True, exist_ok=True)
edges.parent.mkdir(parents=True, exist_ok=True)
current = 7.1738981169369e-25 if reltol >= 1e-2 else 1.0e-33
cathode = current if reltol >= 1e-2 else -current

with terminal.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "point_index", "bias_V", "bias_contact", "contact",
        "current_electron", "current_hole", "electron_minus_hole",
        "current_total", "electron_plus_hole",
        "current_electron_A_per_um", "current_hole_A_per_um",
        "electron_minus_hole_A_per_um", "current_total_A_per_um",
        "electron_plus_hole_A_per_um", "converged", "solver_method",
        "gummel_iterations", "newton_iterations", "handoff_stage",
    ])
    writer.writerow([0, 0, "Anode", "Anode", 0, -current, current, current, -current,
                     0, -current, current, current, -current, 1, "gummel_newton", 0, 8, "newton"])
    writer.writerow([0, 0, "Anode", "Cathode", 0, -cathode, cathode, cathode, -cathode,
                     0, -cathode, cathode, cathode, -cathode, 1, "gummel_newton", 0, 8, "newton"])

with edges.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["current_contact", "edge_id", "node0", "node1", "current_total_A_per_um", "current_electron", "current_hole"])
    writer.writerow(["Anode", 0, 0, 1, current, 0, -current])
    writer.writerow(["Cathode", 1, 0, 1, cathode, 0, -cathode])

vtk_prefix = Path(cfg.get("sweep", {}).get("vtk_prefix", terminal.parent / "probe"))
vtk = vtk_prefix.with_name(vtk_prefix.name + "_0000.vtk")
vtk.write_text(\"\"\"# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0 0 0
1 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0
\"\"\")
""".lstrip()
        )

    def _write_sentaurus_0v_reference_files(self, reference: Path) -> None:
        source = reference / "source"
        source.mkdir(parents=True)
        lines = ["\n"] * 555
        lines[515] = "contact        voltage     electron current    hole current  conduction current\n"
        lines[516] = " Anode        0.000E+00       6.741E-26        -2.936E-25       -2.262E-25\n"
        lines[517] = " Cathode      0.000E+00      -6.741E-26         2.936E-25        2.262E-25\n"
        lines[546] = "contact        voltage     electron current    hole current  conduction current\n"
        lines[547] = " Anode        0.000E+00       8.309E-26        -8.005E-25       -7.174E-25\n"
        lines[548] = " Cathode      0.000E+00      -8.309E-26         8.005E-25        7.174E-25\n"
        (source / "pn2d_0v.log_des.log").write_text("".join(lines))
        (source / "pn2d_0v.plt").write_text(
            """
DF-ISE text

Info {
  version   = 1.0
  type      = xyplot
  datasets  = [
    "time"
    "Cathode OuterVoltage" "Cathode InnerVoltage" "Cathode QuasiFermiPotential" "Cathode DisplacementCurrent" "Cathode eCurrent"
    "Cathode hCurrent" "Cathode TotalCurrent" "Cathode Charge" "Anode OuterVoltage" "Anode InnerVoltage"
    "Anode QuasiFermiPotential" "Anode DisplacementCurrent" "Anode eCurrent" "Anode hCurrent" "Anode TotalCurrent"
    "Anode Charge" ]
}

Data {
      0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00  -6.74115722458708E-26
      2.93597182315862E-25   2.26185610069991E-25   5.28941943767407E-16   0.00000000000000E+00   0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   6.74115645827379E-26  -2.93597174652734E-25  -2.26185610069996E-25
     -5.28941943767407E-16
      0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00  -8.30875544197565E-26
      8.00477366113444E-25   7.17389811693687E-25   5.28941933930105E-16   0.00000000000000E+00   0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   8.30875493000826E-26  -8.00477360993773E-25  -7.17389811693691E-25
     -5.28941933930104E-16
}
""".lstrip()
        )

    def _write_terminal_balance_csv(self, path: Path, rows: list[list[object]]) -> None:
        self._write_csv(path, [
            "contact",
            "current_electron_A_per_um",
            "current_hole_A_per_um",
            "electron_minus_hole_A_per_um",
            "current_total_A_per_um",
            "electron_plus_hole_A_per_um",
            "converged",
            "solver_method",
            "newton_iterations",
            "handoff_stage",
        ], [
            [*row, 1, "gummel_newton", 3, "newton"] for row in rows
        ])

    def _write_contact_edge_csv(self, path: Path, rows: list[list[object]]) -> None:
        self._write_csv(path, [
            "current_contact",
            "edge_id",
            "node0",
            "node1",
            "current_total_A_per_um",
            "psi0",
            "psi1",
            "phin0",
            "phin1",
            "phip0",
            "phip1",
            "electron_continuity_flux",
            "hole_continuity_flux",
            "current_electron",
            "current_hole",
        ], [
            [contact, edge_id, 0, 1, current, 0, 0, 0, 0, 0, 0, 0, 0, current, 0]
            for contact, edge_id, current in rows
        ])

    def _write_csv(self, path: Path, header: list[str], rows: list[list[object]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    def _write_geometry_fixture(
        self,
        export: Path,
        tdx: Path,
        *,
        tdr_triangle: list[int] | None = None,
        tdr_node1_x: float = 1.0,
    ) -> None:
        tdx.mkdir(parents=True)
        self._write_csv(export / "nodes.csv", ["id", "x_um", "y_um"], [
            [0, 0.0, 0.0],
            [1, tdr_node1_x, 0.0],
            [2, 0.0, 1.0],
        ])
        triangle = tdr_triangle or [0, 1, 2]
        self._write_csv(export / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
            [0, *triangle, "R.Si", "Si"],
        ])
        self._write_csv(export / "contacts.csv", ["name", "node_ids", "region"], [
            ["Anode", "0;2", "R.Si"],
        ])
        (export / "metadata.json").write_text(json.dumps({
            "vertex_count": 3,
            "regions": [
                {"name": "R.Si", "material": "Silicon", "triangles": 1, "edges": 0},
                {"name": "Anode", "material": "Contact", "triangles": 0, "edges": 1},
            ],
        }, indent=2) + "\n")
        (tdx / "pn2d_msh.grd").write_text(
            """
DF-ISE text
Info {
  type = grid
  dimension = 2
  nb_vertices = 3
  nb_edges = 3
  nb_elements = 2
  nb_regions = 2
  regions = [ "R.Si" "Anode" ]
  materials = [ Silicon Contact ]
}
Data {
  Vertices (3) {
    0 0
    1 0
    0 1
  }
  Edges (3) {
    0 1
    1 2
    2 0
  }
  Elements (2) {
    2 0 1 2
    1 0 2
  }
  Region ("R.Si") {
    material = Silicon
    Elements (1) { 0 }
  }
  Region ("Anode") {
    material = Contact
    Elements (1) { 1 }
  }
}
""".strip()
            + "\n"
        )


if __name__ == "__main__":
    unittest.main()
