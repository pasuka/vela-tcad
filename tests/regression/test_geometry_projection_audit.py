"""Regression coverage for geometry projection audit diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_geometry_projection.py"


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def compare_row(bias: float, node_id: int, quantity: str, sentaurus: float) -> list[object]:
    return [bias, bias, node_id, float(node_id), 0.0, quantity, "", sentaurus, "", 0.0, "", "", "", ""]


class GeometryProjectionAuditTest(unittest.TestCase):
    def test_reports_unavailable_box_face_normal_and_edge_projection_ratio(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_geometry_projection_") as td:
            tmp = Path(td)
            internal = tmp / "internal.csv"
            edges = tmp / "sg_edges.csv"
            elements = tmp / "elements.csv"
            compare = tmp / "compare.csv"
            cell = tmp / "cell.csv"
            out_csv = tmp / "geometry_projection_audit.csv"
            out_md = tmp / "geometry_projection_audit_summary.md"

            bias = -20.0
            write_csv(
                internal,
                [
                    "source_location_type",
                    "source_entity_id",
                    "bias_V",
                    "Jn_mag_used_A_per_cm2",
                    "Jp_mag_used_A_per_cm2",
                    "source_weight_or_volume_cm2_for_2D",
                ],
                [["edge", 0, bias, 3.0, 4.0, 1.0e-8]],
            )
            write_csv(
                edges,
                [
                    "point_index",
                    "bias_V",
                    "edge_id",
                    "node0",
                    "node1",
                    "x0_um",
                    "y0_um",
                    "x1_um",
                    "y1_um",
                    "edge_length_m",
                    "edge_area_proxy_m2",
                    "electron_raw_signed_flux_proxy",
                    "hole_raw_signed_flux_proxy",
                ],
                [[0, bias, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0e-6, 1.0e-12, -1.0, 1.0]],
            )
            write_csv(elements, ["id", "node0", "node1", "node2"], [[7, 0, 1, 2]])
            header = [
                "bias_V",
                "vela_bias_V",
                "node_id",
                "x_um",
                "y_um",
                "quantity",
                "sentaurus_field",
                "sentaurus_value",
                "vela_field",
                "vela_value_scaled_to_sentaurus_units",
                "comparison_basis",
                "vela_over_sentaurus",
                "diff",
                "abs_diff",
            ]
            rows = []
            for node in (0, 1):
                rows.extend([
                    compare_row(bias, node, "electron_current_density_x", 5.0),
                    compare_row(bias, node, "electron_current_density_y", 12.0),
                    compare_row(bias, node, "hole_current_density_x", 4.0),
                    compare_row(bias, node, "hole_current_density_y", 0.0),
                ])
            write_csv(compare, header, rows)
            write_csv(
                cell,
                ["bias_V", "mode", "qG_cell_contribution_A_per_um", "qG_sent_cell_contribution_A_per_um"],
                [
                    [bias, "current_cell_vector_ls_from_edge_tangent", 1.0, 2.0],
                    [bias, "current_cell_vector_ls_from_box_normal", 1.0, 2.0],
                ],
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--internal-audit",
                    str(internal),
                    "--self-edge-topology",
                    str(edges),
                    "--elements-csv",
                    str(elements),
                    "--node-compare",
                    str(compare),
                    "--edge-to-cell-reconstruction",
                    str(cell),
                    "--biases",
                    "-20",
                    "--out-csv",
                    str(out_csv),
                    "--out-summary",
                    str(out_md),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with out_csv.open(newline="", encoding="utf-8") as handle:
                rows_out = list(csv.DictReader(handle))
            electron = next(row for row in rows_out if row["carrier"] == "electron")
            self.assertEqual(electron["box_face_normal_x"], "")
            self.assertEqual(electron["projection_type_used_by_diagnostic"], "edge_tangent_projection;box_face_normal_unavailable")
            self.assertAlmostEqual(float(electron["sentaurus_J_vector_mag_A_per_cm2"]), 13.0)
            self.assertAlmostEqual(float(electron["sentaurus_J_edge_projection_A_per_cm2"]), 5.0)
            self.assertAlmostEqual(float(electron["current_scalar_signed_A_per_cm2"]), -3.0)
            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("box-face normal unavailable", summary)
            self.assertIn("previous tangent and box-normal LS modes identical: yes", summary)


if __name__ == "__main__":
    unittest.main()
