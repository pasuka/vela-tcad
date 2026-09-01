import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_g3_sentaurus_equation_balance import (
    analyze_endpoint,
    cosine,
    current_scale_from_sg,
)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


class SentaurusEquationBalanceAuditTest(unittest.TestCase):
    def test_flux_scale_and_endpoint_alignment(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            export = root / "export"
            write_csv(export / "nodes.csv", ["id", "x_um", "y_um"], [
                {"id": 1, "x_um": 2.0, "y_um": 3.0},
                {"id": 2, "x_um": 4.0, "y_um": 5.0},
            ])
            write_csv(
                export / "fields" / "eContinuityRhs_region0.csv",
                ["node_id", "component0"],
                [
                    {"node_id": 1, "component0": 1.0e-15},
                    {"node_id": 2, "component0": -2.0e-15},
                ],
            )
            carrier = root / "carrier.csv"
            write_csv(carrier, [
                "node_id", "x", "y", "electron_residual",
                "electron_flux_abs_sum", "electron_recombination",
            ], [
                {
                    "node_id": 1, "x": 2.0, "y": 3.0,
                    "electron_residual": 3.0,
                    "electron_flux_abs_sum": 6.0,
                    "electron_recombination": 0.0,
                },
                {
                    "node_id": 2, "x": 4.0, "y": 5.0,
                    "electron_residual": 4.0,
                    "electron_flux_abs_sum": 8.0,
                    "electron_recombination": 0.0,
                },
            ])
            sg = root / "sg.csv"
            write_csv(sg, [
                "electron_flux", "electron_particle_line_flux_per_m_s",
            ], [
                {"electron_flux": 2.0, "electron_particle_line_flux_per_m_s": 20.0},
                {"electron_flux": -3.0, "electron_particle_line_flux_per_m_s": -30.0},
            ])
            feedback = root / "feedback.json"
            feedback.write_text(json.dumps({
                "points": [{
                    "sentaurus_terminal_A_per_um": 1.0,
                    "variants": {"SSS": {
                        "current_A_per_um": 1.01,
                        "ratio_to_sentaurus": 1.01,
                    }},
                }],
            }), encoding="utf-8")

            scale = current_scale_from_sg(sg)
            self.assertAlmostEqual(scale["continuity_particle_scale_per_m_s"], 10.0)
            endpoint, rows = analyze_endpoint(
                name="fixture",
                sentaurus_export=export,
                vela_carrier_csv=carrier,
                vela_sg_csv=sg,
                state_feedback_summary=feedback,
                selected_nodes=(1, 2),
            )
            self.assertEqual(len(rows), 2)
            self.assertEqual(endpoint["maximum_coordinate_error_um"], 0.0)
            self.assertAlmostEqual(
                endpoint["sentaurus_native"]["seven_node_l2_A"],
                (5.0 ** 0.5) * 1.0e-15,
            )
            self.assertAlmostEqual(
                endpoint["vela_same_sentaurus_state"]["seven_node_l2_scaled"],
                5.0,
            )
            self.assertAlmostEqual(
                endpoint["vela_same_sentaurus_state"][
                    "seven_node_energy_fraction_of_silicon"
                ],
                1.0,
            )

    def test_cosine(self) -> None:
        self.assertAlmostEqual(cosine([1.0, -2.0], [2.0, -4.0]), 1.0)
        self.assertIsNone(cosine([0.0], [1.0]))


if __name__ == "__main__":
    unittest.main()
