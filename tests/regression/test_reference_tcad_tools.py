#!/usr/bin/env python3
"""Regression coverage for neutral reference TCAD CSV conversion tools."""

from __future__ import annotations

import csv
import json
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
            self.assertEqual(iv["doping"][0]["donors"], 1.0e17)
            self.assertEqual(iv["doping"][1]["acceptors"], 1.0e17)
            self.assertEqual(cv["sweep"]["mode"], "cv_quasistatic")
            self.assertEqual(bv["sweep"]["mode"], "bv_reverse")

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

    @staticmethod
    def _write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
