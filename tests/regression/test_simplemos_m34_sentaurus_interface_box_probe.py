from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import unittest

from scripts.run_simplemos_m34_sentaurus_interface_box_probe import (
    deck,
    infer_input_to_debug_local_permutation,
    parse_debug_block,
    parse_debug_info,
    parse_log_stats,
    raw_coefficient,
    vela_positive_coefficient,
)


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m34_sentaurus_interface_box_probe_evidence.json"
)
REPORT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "sentaurus_interface_box_probe"
    / "m34_sentaurus_interface_box_probe_report.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


DEBUG = """
Info {
dimension = 2
nb_vertices = 4
nb_grd_elements = 3
nb_des_elements = 2
}
Measure { # unit = [um^2]
# grd_elem des_elem elem_type
0 1 2 0.1 0.2 0.3
1 0 2 0.4 0.5 0.6
2 -1 1 # contact or interface
}
Coefficients { # unit = [1]
# grd_elem des_elem elem_type
0 1 2 0.0 0.5 0.5
1 0 2 0.5 0.0 0.5
2 -1 1 # contact or interface
}
"""


LOG = """
Number of grid points is 4.
CVPL_AverageBoxMethod = TRUE
NumberOfEdges = 7
NumberOfGeometricalEdges = 5
NumberOfDoubleEdges = 0
NumberOfVertices = 6
NumberOfElements = 2
NumberOfObtuseElements = 0
NumberOfNonDelaunayElements = 0
"""


class SimpleMosM34ProbeTests(unittest.TestCase):
    def test_deck_requests_direct_average_box_debug(self) -> None:
        text = deck()
        self.assertIn("AverageBoxMethod", text)
        self.assertIn("BoxMeasureFromFile(GrdNumbering)", text)
        self.assertIn("Coupled(Iterations=1) { Poisson }", text)
        self.assertNotIn("Electron Hole", text)

    def test_debug_parser_preserves_grid_and_device_numbering(self) -> None:
        info = parse_debug_info(DEBUG)
        self.assertEqual(info["nb_vertices"], 4)
        self.assertEqual(info["nb_grd_elements"], 3)
        measure = parse_debug_block(DEBUG, "Measure")
        coefficient = parse_debug_block(DEBUG, "Coefficients")
        self.assertEqual(measure[0]["des"], 1)
        self.assertNotIn(2, measure)
        self.assertEqual(coefficient[1]["values"], [0.5, 0.0, 0.5])

    def test_log_parser_separates_grid_and_internal_counts(self) -> None:
        stats = parse_log_stats(LOG)
        self.assertEqual(stats["input_grid_points"], 4)
        self.assertEqual(stats["internal_vertices"], 6)
        self.assertEqual(stats["internal_edges"] - stats["geometrical_edges"], 2)
        self.assertTrue(stats["average_box_method"])

    def test_local_coefficient_matches_cotangent_and_positive_fallback(self) -> None:
        right = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        self.assertTrue(math.isclose(raw_coefficient(right, 0), 0.0, abs_tol=1e-15))
        self.assertTrue(math.isclose(raw_coefficient(right, 1), 0.5))
        obtuse = [(0.0, 0.0), (1.0, 0.0), (-0.2, 0.1)]
        raw = raw_coefficient(obtuse, 0)
        self.assertLess(raw, 0.0)
        self.assertGreater(vela_positive_coefficient(obtuse, 0), 0.0)

    def test_local_permutation_is_inferred_from_coefficients(self) -> None:
        points = {0: (0.0, 0.0), 1: (2.0, 0.0), 2: (0.2, 1.0)}
        elements = {0: {"nodes": [0, 1, 2]}}
        raw = [raw_coefficient([points[i] for i in range(3)], local)
               for local in range(3)]
        coefficients = {0: {"values": [raw[1], raw[0], raw[2]]}}
        permutation, error = infer_input_to_debug_local_permutation(
            {0}, elements, points, coefficients
        )
        self.assertEqual(permutation, (1, 0, 2))
        self.assertLess(error, 1e-15)

    def test_frozen_oracle_closes_interface_geometry_without_double_nodes(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report = json.loads(REPORT.read_text(encoding="utf-8"))
        summary = evidence["summary"]
        interface = report["si_oxide_interface"]
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(summary["triangle_records"], 2746)
        self.assertEqual(summary["input_to_debug_local_permutation"], [1, 0, 2])
        self.assertEqual(summary["internal_vertex_delta"], 82)
        self.assertEqual(summary["unique_contact_nodes"], 82)
        self.assertEqual(summary["double_edges"], 0)
        self.assertFalse(summary["explicit_si_oxide_double_nodes_supported"])
        self.assertEqual(interface["edge_count"], 20)
        self.assertLess(
            abs(interface["vela_region_over_sentaurus_region_poisson_min"] - 1.0),
            1e-12,
        )
        self.assertLess(
            abs(interface["vela_region_over_sentaurus_region_poisson_max"] - 1.0),
            1e-12,
        )
        for relative, expected in evidence["artifact_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)
        for relative, expected in evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
