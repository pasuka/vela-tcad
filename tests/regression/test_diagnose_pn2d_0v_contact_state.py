#!/usr/bin/env python3
"""Regression coverage for scripts/diagnose_pn2d_0v_contact_state.py."""

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
SCRIPT = REPO / "scripts" / "diagnose_pn2d_0v_contact_state.py"

K_B_OVER_Q = 8.617333262145e-5
TEMPERATURE_K = 300.0
VT = K_B_OVER_Q * TEMPERATURE_K
NI_CM3 = 1.0e10
NA_CM3 = 1.0e17
ND_CM3 = 1.0e17
PSI_P = -VT * math.log(NA_CM3 / NI_CM3)
PSI_N = VT * math.log(ND_CM3 / NI_CM3)
VBI = VT * math.log(NA_CM3 * ND_CM3 / (NI_CM3 * NI_CM3))

# Two-row strip mesh, x = 0..5; bottom nodes 0-5, top nodes 6-11.
# Anode contact: {0, 6} (p-side); Cathode contact: {5, 11} (n-side).
# First rings: Anode -> {1, 7}; Cathode -> {4, 10}.
# Quasi-neutral plateau nodes after exclusion: p-side {2, 8}, n-side {3, 9}.
P_SIDE_NODES = [0, 1, 2, 6, 7, 8]
N_SIDE_NODES = [3, 4, 5, 9, 10, 11]


def node_psi(node: int) -> float:
    return PSI_P if node in P_SIDE_NODES else PSI_N


class DiagnosePn2d0vContactStateTest(unittest.TestCase):
    def _write_csv(self, path: Path, header: list[str], rows: list[list[object]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    def _write_fixture(self, root: Path, vela_psi_offset: float = 0.0) -> tuple[Path, Path]:
        reference = root / "reference"
        fields = reference / "sim_fields" / "0v" / "fields"
        vela = reference / "vela"
        fields.mkdir(parents=True)
        vela.mkdir(parents=True)

        nodes = [{"id": i, "x": float(i % 6), "y": float(i // 6)} for i in range(12)]
        triangles = []
        tri_id = 0
        for i in range(5):
            triangles.append({"id": tri_id, "node_ids": [i, i + 1, i + 6]})
            tri_id += 1
            triangles.append({"id": tri_id, "node_ids": [i + 1, i + 7, i + 6]})
            tri_id += 1
        (vela / "mesh.json").write_text(json.dumps({
            "nodes": nodes,
            "triangles": triangles,
            "contacts": [
                {"name": "Anode", "node_ids": [0, 6]},
                {"name": "Cathode", "node_ids": [5, 11]},
            ],
        }, indent=2) + "\n")
        self._write_csv(reference / "contacts.csv", ["name", "node_ids", "region"], [
            ["Anode", "0;6", "R.Si"],
            ["Cathode", "5;11", "R.Si"],
        ])
        self._write_csv(reference / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
            [node, 0.0, NA_CM3] if node in P_SIDE_NODES else [node, ND_CM3, 0.0]
            for node in range(12)
        ])

        sentaurus_psi = [node_psi(node) for node in range(12)]
        n_cm3 = [NI_CM3 * math.exp(psi / VT) for psi in sentaurus_psi]
        p_cm3 = [NI_CM3 * math.exp(-psi / VT) for psi in sentaurus_psi]
        self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"],
                        [[node, value] for node, value in enumerate(sentaurus_psi)])
        self._write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"],
                        [[node, 0.0] for node in range(12)])
        self._write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"],
                        [[node, 0.0] for node in range(12)])
        self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"],
                        [[node, value] for node, value in enumerate(n_cm3)])
        self._write_csv(fields / "hDensity_region0.csv", ["node_id", "component0"],
                        [[node, value] for node, value in enumerate(p_cm3)])

        vela_psi = [value + vela_psi_offset for value in sentaurus_psi]

        def block(values: list[float]) -> str:
            return "\n".join(f"{value:.17g}" for value in values)

        vtk = root / "vela_state.vtk"
        vtk.write_text(
            "# vtk DataFile Version 2.0\n"
            "mini\n"
            "ASCII\n"
            "DATASET UNSTRUCTURED_GRID\n"
            "POINTS 12 float\n"
            + "\n".join(f"{node['x']} {node['y']} 0" for node in nodes) + "\n"
            "POINT_DATA 12\n"
            "SCALARS Potential float 1\n"
            "LOOKUP_TABLE default\n"
            + block(vela_psi) + "\n"
            "SCALARS ElectronQuasiFermi float 1\n"
            "LOOKUP_TABLE default\n"
            + block([0.0] * 12) + "\n"
            "SCALARS HoleQuasiFermi float 1\n"
            "LOOKUP_TABLE default\n"
            + block([0.0] * 12) + "\n"
            "SCALARS Electrons float 1\n"
            "LOOKUP_TABLE default\n"
            + block([value * 1.0e6 for value in n_cm3]) + "\n"
            "SCALARS Holes float 1\n"
            "LOOKUP_TABLE default\n"
            + block([value * 1.0e6 for value in p_cm3]) + "\n"
        )
        return reference, vtk

    def _run_script(self, reference: Path, vtk: Path, reports: Path) -> dict:
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--reference-root",
                str(reference),
                "--existing-vtk",
                str(vtk),
                "--output-dir",
                str(reports),
                "--ni-cm3",
                str(NI_CM3),
                "--bgn",
                "none",
            ],
            check=True,
            cwd=REPO,
        )
        self.assertTrue((reports / "pn2d_0v_contact_state.md").is_file())
        return json.loads((reports / "pn2d_0v_contact_state.json").read_text())

    def test_contact_and_first_ring_delta_statistics(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_contact_state_stats_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_fixture(root)
            report = self._run_script(reference, vtk, root / "reports")

            self.assertEqual(report["status"], "pass")
            mesh = report["mesh"]
            self.assertEqual(mesh["contact_node_counts"], {"Anode": 2, "Cathode": 2})
            self.assertEqual(mesh["first_ring_node_counts"], {"Anode": 2, "Cathode": 2})
            self.assertEqual(mesh["contacts_csv_match"], {"Anode": True, "Cathode": True})

            anode = report["contacts"]["Anode"]
            self.assertEqual(anode["sentaurus"]["contact_psi_V"]["points"], 2)
            self.assertAlmostEqual(anode["sentaurus"]["contact_psi_V"]["median"], PSI_P, places=9)
            self.assertEqual(anode["sentaurus"]["first_ring_psi_V"]["points"], 2)
            self.assertAlmostEqual(anode["sentaurus"]["first_ring_psi_V"]["median"], PSI_P, places=9)
            self.assertAlmostEqual(anode["sentaurus"]["contact_to_first_ring_delta_V"], 0.0, places=9)

            cathode = report["contacts"]["Cathode"]
            self.assertAlmostEqual(cathode["vela"]["contact_psi_V"]["median"], PSI_N, places=9)
            self.assertAlmostEqual(cathode["vela"]["contact_to_first_ring_delta_V"], 0.0, places=9)

            for contact in ("Anode", "Cathode"):
                delta = report["contacts"][contact]["delta_vela_minus_sentaurus"]
                self.assertEqual(delta["contact_psi_V"]["points"], 2)
                self.assertEqual(delta["first_ring_psi_V"]["points"], 2)
                self.assertAlmostEqual(delta["contact_psi_V"]["median"], 0.0, places=9)
                self.assertAlmostEqual(delta["first_ring_psi_V"]["median"], 0.0, places=9)

    def test_built_in_potential_estimate_written_to_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_contact_state_vbi_") as tmp:
            root = Path(tmp)
            reference, vtk = self._write_fixture(root)
            report = self._run_script(reference, vtk, root / "reports")

            built_in = report["built_in_potential"]
            self.assertAlmostEqual(built_in["estimate_V"], VBI, places=9)
            self.assertAlmostEqual(built_in["na_cm3"], NA_CM3)
            self.assertAlmostEqual(built_in["nd_cm3"], ND_CM3)
            self.assertEqual(built_in["plateau_node_counts"], {"p_side": 2, "n_side": 2})
            self.assertAlmostEqual(built_in["sentaurus_plateau_delta_V"], VBI, places=9)
            self.assertAlmostEqual(built_in["vela_plateau_delta_V"], VBI, places=9)
            self.assertEqual(report["classification"], "consistent")
            self.assertFalse(any(report["root_cause_flags"].values()))

    def test_contact_psi_offset_classified_as_pinning_mismatch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_contact_state_pinning_") as tmp:
            root = Path(tmp)
            offset = 0.15
            reference, vtk = self._write_fixture(root, vela_psi_offset=offset)
            report = self._run_script(reference, vtk, root / "reports")

            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["classification"], "contact_psi_pinning_mismatch")
            self.assertTrue(report["root_cause_flags"]["contact_psi_pinning_mismatch"])
            self.assertFalse(report["root_cause_flags"]["built_in_potential_mismatch"])
            self.assertFalse(report["root_cause_flags"]["bulk_poisson_state_mismatch"])
            for contact in ("Anode", "Cathode"):
                delta = report["contacts"][contact]["delta_vela_minus_sentaurus"]
                self.assertAlmostEqual(delta["contact_psi_V"]["median"], offset, places=9)
                self.assertAlmostEqual(delta["first_ring_psi_V"]["median"], offset, places=9)
            parity = report["potential_parity"]
            self.assertAlmostEqual(parity["rms_V"], offset, places=9)
            self.assertLess(parity["rms_after_side_shift_V"], 1.0e-9)
            self.assertTrue(report["contact_offset_explains_potential_rms"])


if __name__ == "__main__":
    unittest.main()
