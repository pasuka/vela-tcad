#!/usr/bin/env python3
"""Regression coverage for Sentaurus text import helpers."""

from __future__ import annotations

import importlib.util
import csv
import json
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

    def test_reference_import_config_generates_iv_bv_tree_and_reports(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_reference_") as tmp:
            root = Path(tmp)
            source = root / "pn2d"
            output = root / "out"
            source.mkdir()
            for name in ["pn2d_des.tdr", "pn2d_bv_des.tdr", "pn2d_msh.tdr"]:
                (source / name).write_text("fake tdr placeholder\n")
            (source / "pn2d_sde.cmd").write_text("(define L 2.0)\n(define H 0.5)\n(define Xj 1.0)\n")
            (source / "pn2d_sdevice.cmd").write_text(
                """
File { Grid="pn2d_msh.tdr" Plot="pn2d_des.tdr" Current="pn2d_iv.plt" }
Electrode { { Name="Anode" Voltage=0 } { Name="Cathode" Voltage=0 } }
Physics { Fermi Mobility(DopingDep HighFieldSaturation) Recombination(SRH Auger) EffectiveIntrinsicDensity(OldSlotboom) }
Solve {
  Quasistationary( Goal { Name="Anode" Voltage=1.0 } ){ Coupled { Poisson Electron Hole } }
}
""".strip()
                + "\n"
            )
            (source / "pn2d_bv_sdevice.cmd").write_text(
                """
File { Grid="pn2d_msh.tdr" Plot="pn2d_bv_des.tdr" Current="pn2d_bv.plt" }
Electrode { { Name="Anode" Voltage=0 } { Name="Cathode" Voltage=0 } }
Physics { Fermi Avalanche(OkutoCrowell) }
Solve {
  Quasistationary( Goal { Name="Cathode" Voltage=50.0 } ){ Coupled { Poisson Electron Hole } }
}
""".strip()
                + "\n"
            )
            (source / "pn2d_iv.plt").write_text(
                """
Info { datasets = ["time" "Anode OuterVoltage" "Anode TotalCurrent"] }
Data {
  0 0 0
  1 0.5 1e-9
  2 1.0 1e-6
}
""".strip()
                + "\n"
            )
            (source / "pn2d_bv.plt").write_text(
                """
Info { datasets = ["time" "Cathode OuterVoltage" "Cathode TotalCurrent"] }
Data {
  0 0 0
  1 25 1e-15
  2 50 2e-15
}
""".strip()
                + "\n"
            )
            config = root / "pn2d_reference.json"
            config.write_text(json.dumps({
                "case": "pn2d",
                "device": "pn_diode",
                "mesh_tdr": "pn2d_msh.tdr",
                "sde_cmd": "pn2d_sde.cmd",
                "tdr_doping": {
                    "compensated_node_policy": "dominant_signed_region",
                },
                "vela_solver": {
                    "method": "gummel_newton",
                    "max_iter": 40,
                    "reltol": 1.0e-8,
                    "abstol": 1.0e-18,
                    "damping_psi": 0.25,
                    "damping_factor": 1.0,
                    "line_search": True,
                    "warm_start": True,
                    "contact_boundary_reconstruction": "dominant_signed_contact_mean",
                    "verbose": False,
                    "handoff": {
                        "fallback": "none",
                        "require_gummel_convergence": True,
                    },
                },
                "simulations": [
                    {
                        "name": "iv",
                        "tdr": "pn2d_des.tdr",
                        "cmd": "pn2d_sdevice.cmd",
                        "plt": "pn2d_iv.plt",
                        "bias_column": "Anode OuterVoltage",
                        "current_column": "Anode TotalCurrent",
                        "kind": "iv",
                        "vela_current_contact": "Cathode",
                        "runtime_doping_scale": 1.0e-4,
                        "runtime_step": 0.1,
                        "comparison_min_points": 3,
                        "comparison": {
                            "candidate_column": "current_total_A_per_um",
                        },
                    },
                    {
                        "name": "bv",
                        "tdr": "pn2d_bv_des.tdr",
                        "cmd": "pn2d_bv_sdevice.cmd",
                        "plt": "pn2d_bv.plt",
                        "bias_column": "Cathode OuterVoltage",
                        "current_column": "Cathode TotalCurrent",
                        "kind": "bv",
                        "runtime_doping_scale": 1.0e-4,
                        "runtime_step": 5.0,
                        "comparison_kind": "iv",
                        "require_trend_match": False,
                        "comparison": {
                            "candidate_column": "current_total_A_per_um",
                        },
                    },
                ],
            }) + "\n")
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
    write_csv(out / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0, 0], [1, 1, 0], [2, 2, 0], [3, 0, 0.5], [4, 1, 0.5], [5, 2, 0.5]])
    write_csv(out / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [[0, 0, 1, 4, "R.Si", "Si"], [1, 0, 4, 3, "R.Si", "Si"], [2, 1, 2, 5, "R.NRegion", "Si"], [3, 1, 5, 4, "R.NRegion", "Si"]])
    write_csv(out / "contacts.csv", ["name", "node_ids", "region"], [["Anode", "0;3", "R.Si"], ["Cathode", "2;5", "R.NRegion"]])
    write_csv(out / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [[0, 0, 1e17], [1, 0, 1e17], [2, 1e17, 0], [3, 0, 1e17], [4, 0, 1e17], [5, 1e17, 0]])
    (out / "doping_metadata.json").write_text(json.dumps({"compensated_nodes": {"policy": args.compensated_doping_policy}}, indent=2) + "\\n")
    (out / "metadata.json").write_text(json.dumps({"vertex_count": 6, "region_count": 2, "dataset_count": 0}, indent=2) + "\\n")
    (out / "field_manifest.json").write_text(json.dumps({"fields": []}, indent=2) + "\\n")
if args.inventory_json:
    Path(args.inventory_json).write_text(json.dumps({"vertex_count": 6, "region_count": 2, "dataset_count": 0}, indent=2) + "\\n")
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
out.parent.mkdir(parents=True, exist_ok=True)
kind = cfg["sweep"]["mode"]
with out.open("w", newline="") as handle:
    writer = csv.writer(handle)
    if kind == "bv_reverse":
        writer.writerow(["bias_V", "max_electric_field_V_per_cm", "current_total", "current_total_A_per_um"])
        writer.writerows([[0, 1e4, 0, 0], [25, 2e4, 1e-9, 1e-15], [50, 3e4, 2e-9, 2e-15]])
    else:
        writer.writerow(["bias_V", "current_total", "current_total_A_per_um"])
        writer.writerows([[0, 0, 0], [0.5, 1e-3, 1e-9], [1.0, 1.0, 1e-6]])
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

            for required in [
                "sde_summary.json",
                "reference_tcad_manifest.json",
                "reference_curves/pn2d_iv_reference.csv",
                "reference_curves/pn2d_bv_reference.csv",
                "cmd/iv_summary.json",
                "cmd/bv_summary.json",
                "vela/mesh.json",
                "vela/simulation_iv.json",
                "vela/simulation_iv_runtime.json",
                "vela/simulation_bv.json",
                "vela/simulation_bv_runtime.json",
                "vela/pn2d_iv.csv",
                "vela/pn2d_bv.csv",
                "reports/pn2d_iv_comparison.json",
                "reports/pn2d_bv_comparison.json",
            ]:
                self.assertTrue((output / required).is_file(), required)

            manifest = json.loads((output / "reference_tcad_manifest.json").read_text())
            self.assertEqual(manifest["schema"], "vela.reference_tcad.sentaurus_reference.v1")
            self.assertEqual(manifest["case"], "pn2d")
            self.assertFalse(manifest["commit_policy"]["raw_sentaurus_artifacts"])
            self.assertIn("doping_metadata.json", manifest["generated"])
            doping_metadata = json.loads((output / "doping_metadata.json").read_text())
            self.assertEqual(
                doping_metadata["compensated_nodes"]["policy"],
                "dominant_signed_region",
            )
            self.assertIn("Avalanche", manifest["unsupported_physics"])
            self.assertIn("reports/pn2d_iv_comparison.json", manifest["comparison_reports"])
            self.assertIn("Fermi statistics approximated by Boltzmann carrier statistics", manifest["warnings"])
            self.assertIn(
                "iv runtime deck uses region-average doping scaled by 0.0001 for Vela convergence diagnostics",
                manifest["warnings"],
            )
            self.assertIn(
                "bv runtime deck uses region-average doping scaled by 0.0001 for Vela convergence diagnostics",
                manifest["warnings"],
            )
            self.assertEqual(len(manifest["warnings"]), len(set(manifest["warnings"])))
            iv_deck = json.loads((output / "vela" / "simulation_iv.json").read_text())
            iv_runtime_deck = json.loads((output / "vela" / "simulation_iv_runtime.json").read_text())
            bv_deck = json.loads((output / "vela" / "simulation_bv.json").read_text())
            bv_runtime_deck = json.loads((output / "vela" / "simulation_bv_runtime.json").read_text())
            self.assertEqual(iv_deck["sweep"]["contact"], "Anode")
            self.assertEqual(iv_deck["sweep"]["current_contact"], "Cathode")
            self.assertEqual(iv_deck["sweep"]["stop"], 1.0)
            self.assertEqual(iv_deck["node_doping_file"], "doping.csv")
            self.assertEqual(iv_deck["solver"]["method"], "gummel_newton")
            self.assertEqual(iv_deck["solver"]["max_iter"], 40)
            self.assertEqual(iv_deck["solver"]["reltol"], 1.0e-8)
            self.assertEqual(iv_deck["solver"]["damping_psi"], 0.25)
            self.assertTrue(iv_deck["solver"]["warm_start"])
            self.assertEqual(
                iv_deck["solver"]["contact_boundary_reconstruction"],
                "dominant_signed_contact_mean",
            )
            self.assertEqual(iv_deck["solver"]["handoff"]["fallback"], "none")
            self.assertEqual(iv_deck["solver"]["mobility"]["model"], "caughey_thomas_field")
            self.assertEqual(iv_deck["solver"]["recombination"], ["srh", "auger"])
            self.assertEqual(iv_deck["solver"]["bandgap_narrowing"], "slotboom")
            self.assertNotIn("impact_ionization", iv_deck["solver"])
            self.assertNotIn("runtime_approximation", iv_deck.get("sentaurus_import", {}))
            iv_report = json.loads((output / "reports" / "pn2d_iv_comparison.json").read_text())
            self.assertEqual(iv_report["status"], "pass")
            self.assertEqual(iv_report["iv"]["points_compared"], 3)
            self.assertEqual(iv_report["iv"]["candidate_column"], "current_total_A_per_um")
            self.assertEqual(iv_runtime_deck["solver"]["method"], "gummel")
            self.assertIn("runtime_approximation", iv_runtime_deck["sentaurus_import"])
            self.assertNotIn("node_doping_file", iv_runtime_deck)
            self.assertEqual(iv_runtime_deck["sweep"]["step"], 0.1)
            self.assertEqual(iv_runtime_deck["doping"][0]["acceptors"], 1.0e13)
            self.assertEqual(bv_deck["sweep"]["mode"], "bv_reverse")
            self.assertEqual(bv_deck["sweep"]["contact"], "Cathode")
            self.assertEqual(bv_deck["sweep"]["stop"], 50.0)
            self.assertEqual(bv_deck["solver"]["method"], "gummel_newton")
            self.assertEqual(bv_deck["solver"]["impact_ionization"]["model"], "selberherr")
            self.assertIn("OkutoCrowell approximated by Selberherr", manifest["warnings"])
            self.assertEqual(bv_runtime_deck["solver"]["method"], "gummel")
            self.assertNotIn("node_doping_file", bv_runtime_deck)
            self.assertEqual(bv_runtime_deck["sweep"]["mode"], "bv_reverse")
            self.assertEqual(bv_runtime_deck["sweep"]["step"], 5.0)
            bv_report = json.loads((output / "reports" / "pn2d_bv_comparison.json").read_text())
            self.assertEqual(bv_report["status"], "pass")
            self.assertEqual(bv_report["checked_kinds"], ["iv"])
            self.assertEqual(bv_report["iv"]["candidate_column"], "current_total_A_per_um")

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
            self.assertEqual(summary["sources"]["tdr"], str(project / "n20_des.tdr"))
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
            self.assertEqual(zero_deck["solver"]["handoff"]["fallback"], "none")
            self.assertEqual(zero_deck["solver"]["handoff"]["gummel_max_iter"], 0)
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

    def _write_csv(self, path: Path, header: list[str], rows: list[list[object]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
