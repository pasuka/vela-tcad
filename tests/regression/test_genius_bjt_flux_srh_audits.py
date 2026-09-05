#!/usr/bin/env python3
"""Regression tests for Genius BJT conservative-flux and SRH audits."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def load_script(name: str):
    path = REPO / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FLUX = load_script("audit_genius_bjt_conservative_sections")
SRH = load_script("audit_genius_bjt_srh_alignment")


class GeniusBjtFluxSrhAuditTest(unittest.TestCase):
    def test_section_current_is_orientation_invariant(self) -> None:
        forward = {
            "x0": "0", "x1": "0", "y0": "0", "y1": "1e-6",
            "electron_particle_line_flux_per_m_s": "2",
            "hole_particle_line_flux_per_m_s": "3",
        }
        reverse = {
            "x0": "0", "x1": "0", "y0": "1e-6", "y1": "0",
            "electron_particle_line_flux_per_m_s": "-2",
            "hole_particle_line_flux_per_m_s": "-3",
        }
        for row in (forward, reverse):
            result = FLUX.section_current([row], axis="y", cut_um=0.5)
            self.assertEqual(1, result["crossing_edge_count"])
            self.assertAlmostEqual(-2.0 * FLUX.Q_C * 1.0e-6, result["electron_A_per_um"])
            self.assertAlmostEqual(3.0 * FLUX.Q_C * 1.0e-6, result["hole_A_per_um"])

    def test_section_current_ignores_edges_on_one_side(self) -> None:
        row = {
            "x0": "0", "x1": "0", "y0": "0", "y1": "0.2e-6",
            "electron_particle_line_flux_per_m_s": "2",
            "hole_particle_line_flux_per_m_s": "3",
        }
        result = FLUX.section_current([row], axis="y", cut_um=0.5)
        self.assertEqual(0, result["crossing_edge_count"])
        self.assertEqual(0.0, result["total_A_per_um"])

    def test_parse_scharfetter_parameter_pair(self) -> None:
        text = """Scharfetter * relation\n{\n taumin = 0, 0\n taumax = 1e-5, 3e-6\n Nref = 1e16, 1e16\n gamma = 1, 1\n Talpha = -1.5, -1.5\n Etrap = 0\n}\n"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "models.par"
            path.write_text(text, encoding="utf-8")
            parsed = SRH.parse_scharfetter(path)
        self.assertEqual([1.0e-5, 3.0e-6], parsed["tau_max_s"])
        self.assertEqual([1.0e16, 1.0e16], parsed["reference_doping_cm3"])
        self.assertEqual(0.0, parsed["trap_energy_eV"])

    def test_independent_lumped_area_sums_to_triangle_area(self) -> None:
        mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 2.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
            ],
            "triangles": [{"id": 0, "node_ids": [0, 1, 2]}],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mesh.json"
            path.write_text(json.dumps(mesh), encoding="utf-8")
            areas = SRH.independent_lumped_areas(path)
        self.assertEqual([1.0 / 3.0] * 3, areas)
        self.assertAlmostEqual(1.0, sum(areas))


if __name__ == "__main__":
    unittest.main()
