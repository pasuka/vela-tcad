#!/usr/bin/env python3
"""Regression coverage for neutral reference TCAD CSV conversion tools."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class ReferenceTcadToolsTest(unittest.TestCase):
    def test_pn_export_converts_to_unit_scaling_deck(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_tcad_") as tmp:
            root = Path(tmp)
            export_dir = root / "export"
            out_dir = root / "vela"
            export_dir.mkdir()

            self._write_csv(
                export_dir / "nodes.csv",
                ["id", "x_um", "y_um"],
                [
                    [0, 0.0, 0.0],
                    [1, 1.0, 0.0],
                    [2, 1.0, 1.0],
                    [3, 0.0, 1.0],
                ],
            )
            self._write_csv(
                export_dir / "elements.csv",
                ["id", "node0", "node1", "node2", "region", "material"],
                [
                    [0, 0, 1, 2, "n_region", "Si"],
                    [1, 0, 2, 3, "p_region", "Si"],
                ],
            )
            self._write_csv(
                export_dir / "contacts.csv",
                ["name", "node_ids", "region"],
                [
                    ["anode", "0;3", "p_region"],
                    ["cathode", "1;2", "n_region"],
                ],
            )
            self._write_csv(
                export_dir / "doping.csv",
                ["node_id", "donors_cm3", "acceptors_cm3"],
                [
                    [0, 0.0, 1.0e17],
                    [1, 1.0e17, 0.0],
                    [2, 1.0e17, 0.0],
                    [3, 0.0, 1.0e17],
                ],
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "convert_tcad_export.py"),
                    "--input-dir",
                    str(export_dir),
                    "--output-dir",
                    str(out_dir),
                    "--device",
                    "pn_diode",
                    "--simulation-types",
                    "iv,cv,bv",
                ],
                check=True,
                cwd=REPO,
            )

            mesh = json.loads((out_dir / "mesh.json").read_text())
            iv = json.loads((out_dir / "simulation_iv.json").read_text())
            cv = json.loads((out_dir / "simulation_cv.json").read_text())
            bv = json.loads((out_dir / "simulation_bv.json").read_text())

            self.assertEqual(mesh["nodes"][1]["x"], 1.0)
            self.assertEqual(mesh["triangles"][0]["node_ids"], [0, 1, 2])
            self.assertEqual(mesh["contacts"][0]["node_ids"], [0, 3])
            self.assertEqual(iv["scaling"], {"mode": "unit_scaling"})
            self.assertEqual(iv["mesh_file"], "mesh.json")
            self.assertEqual(iv["node_doping_file"], "doping.csv")
            self.assertEqual(bv["node_doping_file"], "doping.csv")
            self.assertEqual(iv["doping"][0]["donors"], 1.0e17)
            self.assertEqual(iv["doping"][1]["acceptors"], 1.0e17)
            self.assertEqual(iv["doping"][1]["region"], "p_region")
            self.assertEqual(cv["sweep"]["mode"], "cv_quasistatic")
            self.assertEqual(bv["sweep"]["mode"], "bv_reverse")
            self.assertEqual(
                self._read_csv(out_dir / "doping.csv"),
                self._read_csv(export_dir / "doping.csv"),
            )

    def test_compare_reference_curves_writes_json_and_markdown(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_compare_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            report_json = root / "report.json"
            report_md = root / "report.md"

            header = [
                "bias_V",
                "current_total",
                "capacitance_F_per_m",
                "max_electric_field_V_per_cm",
            ]
            self._write_csv(reference, header, [
                [0.0, 1.0e-12, 2.0e-10, 1.0e4],
                [0.5, 1.0e-9, 2.5e-10, 2.0e4],
                [1.0, 1.0e-6, 3.0e-10, 3.0e4],
            ])
            self._write_csv(candidate, header, [
                [0.0, 1.1e-12, 2.1e-10, 1.1e4],
                [0.5, 1.2e-9, 2.6e-10, 2.1e4],
                [1.0, 1.1e-6, 3.1e-10, 3.1e4],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(candidate),
                    "--output-json",
                    str(report_json),
                    "--output-md",
                    str(report_md),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads(report_json.read_text())
            self.assertTrue(report["iv"]["trend_match"])
            self.assertTrue(report["cv"]["trend_match"])
            self.assertTrue(report["bv"]["trend_match"])
            self.assertIn("orders_of_magnitude", report["iv"])
            self.assertIn("Reference TCAD Curve Comparison", report_md.read_text())

    def test_compare_reference_curves_enforces_single_curve_thresholds(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_compare_gate_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            bad_candidate = root / "bad_candidate.csv"
            report_json = root / "report.json"
            report_md = root / "report.md"

            header = ["bias_V", "current_total"]
            self._write_csv(reference, header, [
                [0.0, 1.0e-12],
                [0.5, 1.0e-9],
                [1.0, 1.0e-6],
            ])
            self._write_csv(candidate, header, [
                [0.0, 1.0e-12],
                [0.5, 1.1e-9],
                [1.0, 1.05e-6],
            ])
            self._write_csv(bad_candidate, header, [
                [0.0, 1.0e-6],
                [0.5, 1.0e-9],
                [1.0, 1.0e-12],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(candidate),
                    "--output-json",
                    str(root / "default_report.json"),
                    "--output-md",
                    str(root / "default_report.md"),
                ],
                check=True,
                cwd=REPO,
            )
            default_report = json.loads((root / "default_report.json").read_text())
            self.assertEqual(default_report["status"], "pass")
            self.assertFalse(default_report["cv"]["available"])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(candidate),
                    "--output-json",
                    str(report_json),
                    "--output-md",
                    str(report_md),
                    "--kind",
                    "iv",
                    "--require-trend-match",
                    "--min-points",
                    "3",
                    "--max-relative-error",
                    "0.2",
                    "--max-orders-of-magnitude",
                    "0.1",
                ],
                check=True,
                cwd=REPO,
            )
            report = json.loads(report_json.read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["iv"]["points_compared"], 3)
            self.assertEqual(report["checked_kinds"], ["iv"])

            failed = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(bad_candidate),
                    "--output-json",
                    str(root / "bad_report.json"),
                    "--output-md",
                    str(root / "bad_report.md"),
                    "--kind",
                    "iv",
                    "--require-trend-match",
                    "--min-points",
                    "3",
                    "--max-orders-of-magnitude",
                    "0.1",
                ],
                cwd=REPO,
            )
            self.assertNotEqual(failed.returncode, 0)
            bad_report = json.loads((root / "bad_report.json").read_text())
            self.assertEqual(bad_report["status"], "fail")
            self.assertTrue(any("trend" in failure for failure in bad_report["failures"]))

    def test_checked_in_pn_validation_assets_are_complete(self) -> None:
        pn_dir = REPO / "reference_tcad" / "pn_diode"
        vela_dir = pn_dir / "vela"
        reports_dir = pn_dir / "reports"
        validation_doc = REPO / "docs" / "validation" / "pn_diode_unit_scaling_validation.md"

        readme = (pn_dir / "README.md").read_text()
        for required in [
            "2D silicon",
            "p_region",
            "n_region",
            "anode",
            "cathode",
            "abrupt junction",
            "forward IV",
            "reverse quasi-static CV",
            "reverse BV",
            "unit_scaling",
        ]:
            self.assertIn(required, readme)

        for csv_name in ["nodes.csv", "elements.csv", "contacts.csv", "doping.csv"]:
            self.assertTrue((pn_dir / csv_name).is_file(), csv_name)

        mesh = json.loads((vela_dir / "mesh.json").read_text())
        self.assertEqual({region["name"] for region in mesh["regions"]}, {"n_region", "p_region"})
        self.assertEqual({contact["name"] for contact in mesh["contacts"]}, {"anode", "cathode"})

        for deck_name in ["simulation_iv.json", "simulation_cv.json", "simulation_bv.json"]:
            deck = json.loads((vela_dir / deck_name).read_text())
            self.assertEqual(deck["scaling"], {"mode": "unit_scaling"})
            self.assertEqual(deck["mesh_file"], "mesh.json")

        cv_deck = json.loads((vela_dir / "simulation_cv.json").read_text())
        bv_deck = json.loads((vela_dir / "simulation_bv.json").read_text())
        self.assertLess(cv_deck["sweep"]["stop"], 0.0)
        self.assertLess(bv_deck["sweep"]["stop"], 0.0)

        report = json.loads((reports_dir / "pn_diode_comparison.json").read_text())
        for section in ["iv", "cv", "bv"]:
            self.assertTrue(report[section]["available"], section)
            self.assertTrue(report[section]["trend_match"], section)

        doc = validation_doc.read_text()
        for required in [
            "IV monotonic",
            "finite capacitance",
            "max field non-decreasing",
            "trend and order-of-magnitude",
            "no calibration claim",
        ]:
            self.assertIn(required, doc)

    def test_checked_in_mos_validation_assets_are_complete(self) -> None:
        cases = [
            {
                "device": "nmos2d",
                "regions": {"p_body", "n_source", "n_drain", "gate_oxide"},
                "idvg_stop_sign": 1.0,
                "idvd_current_sign": 1.0,
                "validation_doc": "nmos2d_unit_scaling_validation.md",
            },
            {
                "device": "pmos2d",
                "regions": {"n_body", "p_source", "p_drain", "gate_oxide"},
                "idvg_stop_sign": -1.0,
                "idvd_current_sign": 1.0,
                "validation_doc": "pmos2d_unit_scaling_validation.md",
            },
        ]

        for case in cases:
            with self.subTest(device=case["device"]):
                device_dir = REPO / "reference_tcad" / case["device"]
                vela_dir = device_dir / "vela"
                reports_dir = device_dir / "reports"

                readme = (device_dir / "README.md").read_text()
                for required in [
                    "mixed Si/SiO2",
                    "metal_gate",
                    "interface charge",
                    "surface mobility",
                    "multi-terminal quasi-static CV",
                    "all-Si MOS baseline",
                    "unit_scaling",
                ]:
                    self.assertIn(required, readme)

                for csv_name in ["nodes.csv", "elements.csv", "contacts.csv", "doping.csv"]:
                    self.assertTrue((device_dir / csv_name).is_file(), csv_name)

                mesh = json.loads((vela_dir / "mesh.json").read_text())
                self.assertEqual({region["name"] for region in mesh["regions"]}, case["regions"])
                self.assertEqual(
                    {contact["name"] for contact in mesh["contacts"]},
                    {"source", "drain", "body", "gate"},
                )
                self.assertIn("SiO2", {region["material"] for region in mesh["regions"]})

                for deck_name in [
                    "simulation_idvd.json",
                    "simulation_idvg.json",
                    "simulation_idvg_surface.json",
                    "simulation_cv.json",
                    "simulation_bv.json",
                ]:
                    deck = json.loads((vela_dir / deck_name).read_text())
                    self.assertEqual(deck["scaling"], {"mode": "unit_scaling"})
                    self.assertEqual(deck["mesh_file"], "mesh.json")

                surface_deck = json.loads((vela_dir / "simulation_idvg_surface.json").read_text())
                self.assertEqual(surface_deck["solver"]["mobility"]["model"], "caughey_thomas_field_surface")
                self.assertEqual(surface_deck["contacts"][2]["type"], "metal_gate")
                self.assertIn("interfaces", surface_deck)

                idvg_rows = self._read_csv(vela_dir / f"{case['device']}_idvg.csv")
                idvg_currents = [abs(float(row["current_total"])) for row in idvg_rows]
                self.assertGreater(idvg_currents[-1], idvg_currents[0])
                self.assertEqual(
                    self._sign(float(idvg_rows[-1]["bias_V"])),
                    case["idvg_stop_sign"],
                )

                idvd_rows = self._read_csv(vela_dir / f"{case['device']}_idvd.csv")
                nonzero_idvd = [float(row["current_total"]) for row in idvd_rows if abs(float(row["bias_V"])) > 0.0]
                self.assertTrue(nonzero_idvd)
                self.assertEqual(self._sign(nonzero_idvd[-1]), case["idvd_current_sign"])

                cv_rows = self._read_csv(vela_dir / f"{case['device']}_cv.csv")
                cv_columns = set(cv_rows[0])
                for column in [
                    "charge_gate_C_per_m",
                    "capacitance_Cgate_gate_F_per_m",
                    "charge_drain_C_per_m",
                    "capacitance_Cgate_drain_F_per_m",
                    "charge_source_C_per_m",
                    "capacitance_Cgate_source_F_per_m",
                    "charge_body_C_per_m",
                    "capacitance_Cgate_body_F_per_m",
                ]:
                    self.assertIn(column, cv_columns)
                    values = [float(row[column]) for row in cv_rows]
                    self.assertTrue(all(math.isfinite(value) for value in values), column)

                bv_rows = self._read_csv(vela_dir / f"{case['device']}_bv.csv")
                fields = [float(row["max_electric_field_V_per_cm"]) for row in bv_rows]
                self.assertTrue(all(b >= a for a, b in zip(fields, fields[1:])))

                for report_name in ["idvd", "idvg", "cv", "bv"]:
                    report = json.loads((reports_dir / f"{case['device']}_{report_name}_comparison.json").read_text())
                    section = "iv" if report_name in {"idvd", "idvg"} else report_name
                    self.assertTrue(report[section]["available"], report_name)
                    self.assertTrue(report[section]["trend_match"], report_name)

                doc = (REPO / "docs" / "validation" / case["validation_doc"]).read_text()
                for required in [
                    "Id-Vg",
                    "Id-Vd",
                    "multi-terminal CV",
                    "BV max field",
                    "trend and order-of-magnitude",
                    "no calibration claim",
                ]:
                    self.assertIn(required, doc)

    def test_checked_in_power_device_validation_assets_are_complete(self) -> None:
        cases = [
            {
                "device": "ldmos2d",
                "regions": {"p_body", "n_source", "n_drift", "n_drain", "gate_oxide"},
                "contacts": {"source", "body", "drain", "gate"},
                "decks": [
                    "simulation_iv.json",
                    "simulation_bv.json",
                    "simulation_bv_fieldplate.json",
                ],
                "reports": ["iv", "bv", "fieldplate"],
                "validation_doc": "ldmos2d_unit_scaling_validation.md",
                "required_doc": [
                    "low-bias DD-IV",
                    "field-plate",
                    "max field",
                    "engineering trend validation",
                    "no calibration claim",
                ],
            },
            {
                "device": "igbt2d",
                "regions": {"p_collector", "n_buffer", "n_drift", "p_base", "n_emitter"},
                "contacts": {"collector", "gate", "emitter"},
                "decks": [
                    "simulation_iv.json",
                    "simulation_high_injection_iv.json",
                    "simulation_charge_cv.json",
                    "simulation_bv.json",
                    "simulation_bv_ii.json",
                ],
                "reports": ["iv", "high_injection_iv", "charge_cv", "bv", "bv_ii"],
                "validation_doc": "igbt2d_unit_scaling_validation.md",
                "required_doc": [
                    "high-injection",
                    "stored charge",
                    "breakdown diagnostic",
                    "engineering trend validation",
                    "no calibration claim",
                ],
            },
        ]

        for case in cases:
            with self.subTest(device=case["device"]):
                device_dir = REPO / "reference_tcad" / case["device"]
                vela_dir = device_dir / "vela"
                reports_dir = device_dir / "reports"
                reference_dir = device_dir / "reference_curves"

                readme = (device_dir / "README.md").read_text()
                for required in [
                    "unit_scaling",
                    "engineering trend validation",
                    "no calibration claim",
                ]:
                    self.assertIn(required, readme)

                for csv_name in ["nodes.csv", "elements.csv", "contacts.csv", "doping.csv"]:
                    self.assertTrue((device_dir / csv_name).is_file(), csv_name)

                mesh = json.loads((vela_dir / "mesh.json").read_text())
                self.assertEqual({region["name"] for region in mesh["regions"]}, case["regions"])
                self.assertEqual({contact["name"] for contact in mesh["contacts"]}, case["contacts"])

                for deck_name in case["decks"]:
                    deck = json.loads((vela_dir / deck_name).read_text())
                    self.assertEqual(deck["scaling"], {"mode": "unit_scaling"})
                    self.assertEqual(deck["mesh_file"], "mesh.json")

                for report_name in case["reports"]:
                    self.assertTrue((reference_dir / f"{case['device']}_{report_name}_reference.csv").is_file())
                    report = json.loads(
                        (reports_dir / f"{case['device']}_{report_name}_comparison.json").read_text()
                    )
                    available = [
                        key for key in ("iv", "cv", "bv")
                        if report[key]["available"] and report[key]["trend_match"]
                    ]
                    self.assertTrue(available, report_name)

                doc = (REPO / "docs" / "validation" / case["validation_doc"]).read_text()
                for required in case["required_doc"]:
                    self.assertIn(required, doc)

    @staticmethod
    def _write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    @staticmethod
    def _read_csv(path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))

    @staticmethod
    def _sign(value: float) -> float:
        if value > 0.0:
            return 1.0
        if value < 0.0:
            return -1.0
        return 0.0


if __name__ == "__main__":
    unittest.main()
