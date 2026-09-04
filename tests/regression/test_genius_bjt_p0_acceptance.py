#!/usr/bin/env python3
"""Regression coverage for the Genius BJT P0 accepted-state contract."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_genius_bjt_accepted_state_pipeline.py"
SPEC = importlib.util.spec_from_file_location("genius_bjt_p0", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
WORKFLOW = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WORKFLOW)


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class GeniusBjtP0AcceptanceTest(unittest.TestCase):
    def test_common_config_enforces_physical_acceptance_gates(self) -> None:
        solver = WORKFLOW.common_config()["solver"]
        rows = solver["carrier_row_convergence"]
        self.assertEqual("enforce", rows["mode"])
        self.assertEqual(1.0e16, rows["min_carrier_density_m3"])
        self.assertEqual(1.0e-6, rows["min_flux_scale_fraction"])
        self.assertEqual(1.0e-6, rows["min_source_global_fraction"])
        self.assertEqual("enforce", solver["global_continuity_closure"]["mode"])

    def test_post_accept_probe_uses_carrier_driving_potential_update(self) -> None:
        state_rows = [
            {"node_id": "0", "electrons_m3": "1e16", "holes_m3": "1e16"},
            {"node_id": "1", "electrons_m3": "1", "holes_m3": "1"},
        ]
        probe_rows = [
            {"node_id": 0, "delta_psi_V": 0.01, "delta_phin_V": 0.04, "delta_phip_V": -0.01},
            {"node_id": 1, "delta_psi_V": 0.0, "delta_phin_V": 10.0, "delta_phip_V": -10.0},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "step.csv"
            write_rows(path, probe_rows)
            summary = WORKFLOW.step_summary(path, state_rows)

        self.assertAlmostEqual(0.03, summary["max_abs_delta_relevant_V"])
        self.assertAlmostEqual(10.0, summary["max_abs_delta_phin_minus_psi_V"])
        self.assertAlmostEqual(10.0, summary["max_abs_delta_phip_minus_psi_V"])


if __name__ == "__main__":
    unittest.main()
