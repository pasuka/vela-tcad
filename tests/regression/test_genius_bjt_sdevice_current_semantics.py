from __future__ import annotations

import importlib.util
import csv
import math
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "diagnose_genius_bjt_sdevice_current_semantics.py"
)
SPEC = importlib.util.spec_from_file_location("genius_bjt_sdevice_current_semantics", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GeniusBjtSdeviceCurrentSemanticsTest(unittest.TestCase):
    def test_area_reconstruction_weights_incident_cells(self) -> None:
        recovered = MODULE.reconstruct(
            4,
            [(0, 1, 2), (0, 2, 3)],
            [1.0, 3.0],
            [(2.0, 4.0), (6.0, 8.0)],
            "area",
        )
        self.assertEqual(recovered[0], (5.0, 7.0))
        self.assertEqual(recovered[1], (2.0, 4.0))
        self.assertEqual(recovered[2], (5.0, 7.0))
        self.assertEqual(recovered[3], (6.0, 8.0))

    def test_least_squares_scale_recovers_global_unit_factor(self) -> None:
        reference = [(1.0, 2.0), (-3.0, 4.0)]
        candidate = [(1.0e6, 2.0e6), (-3.0e6, 4.0e6)]
        scale = MODULE.least_squares_scale(reference, candidate, [0, 1])
        self.assertTrue(math.isclose(scale, 1.0e-6, rel_tol=1.0e-14))

    def test_field_metrics_reports_exact_identity(self) -> None:
        values = [(1.0, 2.0), (-3.0, 4.0)]
        metrics = MODULE.field_metrics(values, values, [0, 1])
        self.assertEqual(metrics["magnitude_abs_dex_max"], 0.0)
        self.assertEqual(metrics["normalized_vector_rmse"], 0.0)
        self.assertEqual(metrics["cosine_similarity"], 1.0)

    def test_sg_line_flux_recovers_cell_current_density(self) -> None:
        charge = 1.602176634e-19
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "edges.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(
                    stream,
                    fieldnames=(
                        "node0", "node1", "x0", "y0", "x1", "y1",
                        "length_m", "couple_m", "hole_particle_line_flux_per_m_s",
                    ),
                )
                writer.writeheader()
                writer.writerows([
                    {
                        "node0": 0, "node1": 1, "x0": 0, "y0": 0,
                        "x1": 1, "y1": 0, "length_m": 1, "couple_m": 2,
                        "hole_particle_line_flux_per_m_s": 2 * 2 / (charge * 1.0e-4),
                    },
                    {
                        "node0": 0, "node1": 2, "x0": 0, "y0": 0,
                        "x1": 0, "y1": 1, "length_m": 1, "couple_m": 3,
                        "hole_particle_line_flux_per_m_s": 3 * 3 / (charge * 1.0e-4),
                    },
                    {
                        "node0": 1, "node1": 2, "x0": 1, "y0": 0,
                        "x1": 0, "y1": 1, "length_m": math.sqrt(2), "couple_m": 0,
                        "hole_particle_line_flux_per_m_s": 0,
                    },
                ])
            values, valid = MODULE.vela_cell_currents_from_sg_edges(path, [(0, 1, 2)])
        self.assertEqual(valid, [0])
        self.assertTrue(math.isclose(values[0][0], 2.0, rel_tol=1.0e-14))
        self.assertTrue(math.isclose(values[0][1], 3.0, rel_tol=1.0e-14))


if __name__ == "__main__":
    unittest.main()
