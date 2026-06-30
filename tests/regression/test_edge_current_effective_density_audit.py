"""Regression coverage for edge current effective-density audit."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_edge_current_effective_density.py"
Q = 1.602176634e-19


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def compare_row(
    bias: float,
    node_id: int,
    x_um: float,
    y_um: float,
    quantity: str,
    sentaurus: float,
    self_value: float,
) -> list[object]:
    return [
        bias,
        bias,
        node_id,
        x_um,
        y_um,
        quantity,
        "",
        sentaurus,
        "",
        self_value,
        "",
        "",
        "",
        "",
    ]


class EdgeCurrentEffectiveDensityAuditTest(unittest.TestCase):
    def test_infers_effective_density_and_best_density_mean(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_edge_eff_density_") as td:
            tmp = Path(td)
            internal = tmp / "internal.csv"
            edges = tmp / "sg_edges.csv"
            compare = tmp / "node_compare.csv"
            out_csv = tmp / "edge_current_effective_density_audit.csv"
            out_md = tmp / "edge_current_effective_density_audit_summary.md"

            bias = -20.0
            length_cm = 1.0e-4
            mu = 1000.0
            phi0 = 0.0
            phi1 = -0.1
            density0 = 1.0e10
            density1 = 9.0e10
            arithmetic = 0.5 * (density0 + density1)
            self_eff = (density0 * density1) ** 0.5
            sent_eff = arithmetic
            self_j = Q * mu * self_eff * abs(phi1 - phi0) / length_cm
            sent_j = Q * mu * sent_eff * abs(phi1 - phi0) / length_cm

            write_csv(
                internal,
                [
                    "source_location_type",
                    "source_entity_id",
                    "bias_V",
                    "x_um",
                    "y_um",
                    "Fn_used_V_per_cm",
                    "Fp_used_V_per_cm",
                    "alpha_n_used_cm_inv",
                    "alpha_p_used_cm_inv",
                    "Jn_mag_used_A_per_cm2",
                    "Jp_mag_used_A_per_cm2",
                ],
                [["edge", 0, bias, 0.5, 0.0, 1000.0, 1000.0, 1.0, 1.0, self_j, self_j]],
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
                ],
                [[0, bias, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0e-6, 1.0e-12]],
            )
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
            rows: list[list[object]] = []
            for node_id, x_um, phi, density in ((0, 0.0, phi0, density0), (1, 1.0, phi1, density1)):
                rows.extend([
                    compare_row(bias, node_id, x_um, 0.0, "electron_qf", phi, phi),
                    compare_row(bias, node_id, x_um, 0.0, "hole_qf", phi, phi),
                    compare_row(bias, node_id, x_um, 0.0, "electron_density", density, density),
                    compare_row(bias, node_id, x_um, 0.0, "hole_density", density, density),
                    compare_row(bias, node_id, x_um, 0.0, "electron_mobility", mu, mu),
                    compare_row(bias, node_id, x_um, 0.0, "hole_mobility", mu, mu),
                    compare_row(bias, node_id, x_um, 0.0, "electron_current_density_x", sent_j, 0.0),
                    compare_row(bias, node_id, x_um, 0.0, "electron_current_density_y", 0.0, 0.0),
                    compare_row(bias, node_id, x_um, 0.0, "hole_current_density_x", sent_j, 0.0),
                    compare_row(bias, node_id, x_um, 0.0, "hole_current_density_y", 0.0, 0.0),
                ])
            write_csv(compare, header, rows)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--internal-audit",
                    str(internal),
                    "--self-edge-topology",
                    str(edges),
                    "--node-compare",
                    str(compare),
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
            self.assertAlmostEqual(float(electron["density_eff_self_from_J"]) / self_eff, 1.0)
            self.assertAlmostEqual(float(electron["density_eff_sent_vector_from_J"]) / sent_eff, 1.0)
            self.assertEqual(electron["best_density_mean_for_self"], "geometric")
            self.assertEqual(electron["best_density_mean_for_sent_vector"], "arithmetic")
            self.assertAlmostEqual(float(electron["self_eff_over_sent_vector_eff"]), self_eff / sent_eff)
            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("self inferred effective density vs Sentaurus", summary)
            self.assertIn("best density mean for Sentaurus vector current", summary)


if __name__ == "__main__":
    unittest.main()
