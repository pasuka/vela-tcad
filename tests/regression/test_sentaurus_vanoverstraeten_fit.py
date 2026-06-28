#!/usr/bin/env python3
"""Regression coverage for the Sentaurus VanOverstraeten fit diagnostic."""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
BUILD_DIR = Path(os.environ.get("VELA_BUILD_DIR", REPO / "build"))
EXE = BUILD_DIR / (
    "sentaurus_vanoverstraeten_fit.exe"
    if sys.platform == "win32"
    else "sentaurus_vanoverstraeten_fit"
)


class SentaurusVanOverstraetenFitTest(unittest.TestCase):
    def test_tool_writes_fit_rows_and_parameter_summary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_vo_fit_") as td:
            tmp = Path(td)
            compare_csv = tmp / "coarse_node_field_compare_aligned.csv"
            elements_csv = tmp / "elements.csv"
            out_csv = tmp / "sentaurus_vanoverstraeten_fit.csv"
            summary_md = tmp / "sentaurus_vanoverstraeten_fit_summary.md"

            elements_csv.write_text(
                "\n".join([
                    "id,node0,node1,node2,region,material",
                    "0,0,1,3,R.Si,Si",
                    "1,1,4,3,R.Si,Si",
                    "2,1,2,4,R.Si,Si",
                    "3,2,5,4,R.Si,Si",
                    "4,3,4,6,R.Si,Si",
                    "5,4,7,6,R.Si,Si",
                    "6,4,5,7,R.Si,Si",
                    "7,5,8,7,R.Si,Si",
                ]),
                encoding="utf-8",
            )

            nodes = [
                (0, 0.0, 0.0),
                (1, 0.5, 0.0),
                (2, 1.0, 0.0),
                (3, 0.0, 0.5),
                (4, 0.5, 0.5),
                (5, 1.0, 0.5),
                (6, 0.0, 1.0),
                (7, 0.5, 1.0),
                (8, 1.0, 1.0),
            ]
            fieldnames = [
                "bias_V",
                "node_id",
                "x_um",
                "y_um",
                "quantity",
                "sentaurus_value",
                "vela_value_scaled_to_sentaurus_units",
            ]

            rows: list[dict[str, str]] = []
            for bias, scale in [("-5", 1.0), ("-10", 1.6)]:
                qf_e: dict[int, float] = {}
                qf_h: dict[int, float] = {}
                for node_id, x_um, y_um in nodes:
                    qf_e[node_id] = scale * (12.0 * x_um * x_um + 8.0 * y_um)
                    qf_h[node_id] = -scale * (15.0 * x_um * x_um + 6.0 * y_um)

                for node_id, x_um, y_um in nodes:
                    f_e = local_field(node_id, nodes, elements_csv, qf_e)
                    f_h = local_field(node_id, nodes, elements_csv, qf_h)
                    alpha_e = 7.0e5 * math.exp(-1.2e6 / max(f_e, 1.0))
                    alpha_h = 1.4e6 * math.exp(-1.8e6 / max(f_h, 1.0))
                    for quantity, sentaurus, self_value in [
                        ("electron_qf", qf_e[node_id], qf_e[node_id]),
                        ("hole_qf", qf_h[node_id], qf_h[node_id]),
                        ("electron_alpha_avalanche", alpha_e, 1.0),
                        ("hole_alpha_avalanche", alpha_h, 1.0),
                    ]:
                        rows.append({
                            "bias_V": bias,
                            "node_id": str(node_id),
                            "x_um": str(x_um),
                            "y_um": str(y_um),
                            "quantity": quantity,
                            "sentaurus_value": f"{sentaurus:.15g}",
                            "vela_value_scaled_to_sentaurus_units": f"{self_value:.15g}",
                        })

            with compare_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)

            result = subprocess.run(
                [
                    str(EXE),
                    "--input-csv",
                    str(compare_csv),
                    "--elements-csv",
                    str(elements_csv),
                    "--output-csv",
                    str(out_csv),
                    "--summary-md",
                    str(summary_md),
                    "--alpha-floor",
                    "1e-30",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with out_csv.open(encoding="utf-8") as handle:
                out_rows = list(csv.DictReader(handle))

            self.assertTrue(out_rows)
            first = out_rows[0]
            for column in [
                "carrier",
                "fit_scope",
                "field_region",
                "A_eff_cm_inv",
                "B_eff_V_per_cm",
                "self_A_cm_inv",
                "self_B_V_per_cm",
                "A_ratio_self_over_fit",
                "B_ratio_self_over_fit",
            ]:
                self.assertIn(column, first)

            scopes = {(row["carrier"], row["fit_scope"], row["field_region"]) for row in out_rows}
            self.assertIn(("electron", "all_biases", "low_mid"), scopes)
            self.assertIn(("hole", "per_bias", "low_mid"), scopes)
            positive_fit = [
                row for row in out_rows
                if row["fit_scope"] == "all_biases" and row["field_region"] == "low_mid"
            ]
            self.assertTrue(any(float(row["A_eff_cm_inv"]) > 0.0 for row in positive_fit))
            self.assertTrue(any(float(row["B_eff_V_per_cm"]) > 0.0 for row in positive_fit))

            summary = summary_md.read_text(encoding="utf-8")
            self.assertIn("F units: V/cm", summary)
            self.assertIn("alpha units: 1/cm", summary)
            self.assertIn("\u8282\u70b9\u8f93\u51fa\u7b49\u6548\u62df\u5408", summary)
            self.assertIn("\u4f4e\u4e2d\u573a", summary)
            self.assertIn("\u9ad8\u573a", summary)


def local_field(
    node_id: int,
    nodes: list[tuple[int, float, float]],
    elements_csv: Path,
    values: dict[int, float],
) -> float:
    node_map = {node[0]: node for node in nodes}
    neighbors: dict[int, set[int]] = {}
    with elements_csv.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            ids = [int(row["node0"]), int(row["node1"]), int(row["node2"])]
            for src in ids:
                neighbors.setdefault(src, set()).update(dst for dst in ids if dst != src)
    _, x0, y0 = node_map[node_id]
    v0 = values[node_id]
    sxx = sxy = syy = sxv = syv = 0.0
    for neighbor in neighbors[node_id]:
        _, x1, y1 = node_map[neighbor]
        dx = x1 - x0
        dy = y1 - y0
        distance = math.hypot(dx, dy)
        if distance <= 0.0:
            continue
        weight = 1.0 / distance
        dv = values[neighbor] - v0
        sxx += weight * dx * dx
        sxy += weight * dx * dy
        syy += weight * dy * dy
        sxv += weight * dx * dv
        syv += weight * dy * dv
    det = sxx * syy - sxy * sxy
    if abs(det) <= 1.0e-60:
        return 0.0
    gx = (sxv * syy - syv * sxy) / det
    gy = (sxx * syv - sxy * sxv) / det
    return math.hypot(gx, gy) * 1.0e4


if __name__ == "__main__":
    unittest.main()
