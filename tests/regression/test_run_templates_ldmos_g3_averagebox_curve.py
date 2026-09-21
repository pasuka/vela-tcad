import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.run_templates_ldmos_g3_averagebox_curve import prepare


class G3AverageBoxCurveTest(unittest.TestCase):
    def test_prepare_freezes_curve_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = {
                "simulation_type": "dc_sweep",
                "solver": {
                    "max_iter": 50,
                    "mobility": {"model": "constant_field"},
                    "impact_ionization": {"model": "none"},
                },
                "sweep": {
                    "bias_points": [float(index) for index in range(31)],
                    "diagnostics": {
                        "terminal_balance": {},
                        "srh_balance": {},
                        "newton_history": {},
                    },
                },
            }
            baseline_path = root / "baseline.json"
            baseline_path.write_text(json.dumps(baseline), encoding="utf-8")
            profile = root / "couples.csv"
            with profile.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["node0", "node1", "couple_m"])
                writer.writerow([0, 1, 2.0e-9])
                writer.writerow([1, 2, 0.0])

            config, count = prepare(baseline_path, profile, root / "run", 400)

            self.assertEqual(count, 2)
            self.assertEqual(len(config["sweep"]["bias_points"]), 31)
            self.assertEqual(config["solver"]["max_iter"], 400)
            self.assertEqual(
                config["solver"]["contact_boundary_reconstruction"],
                "legacy_node_local",
            )
            geometry = config["mesh_geometry"]
            self.assertEqual(geometry["node_volume_policy"], "barycentric")
            self.assertEqual(
                geometry["carrier_transport_couple_profile"],
                "templates_ldmos_external_averagebox",
            )
            self.assertEqual(geometry["external_averagebox_expected_edges"], 2)
            self.assertNotIn("predictor", config["sweep"])
            self.assertEqual(
                config["discretization"]["poisson_charge_volume_policy"],
                "global",
            )

            material_local, _ = prepare(
                baseline_path,
                profile,
                root / "material-local",
                400,
                "material_local",
            )
            self.assertEqual(
                material_local["discretization"][
                    "poisson_charge_volume_policy"
                ],
                "material_local",
            )


if __name__ == "__main__":
    unittest.main()
