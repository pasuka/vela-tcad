import unittest
from pathlib import Path

from scripts.run_templates_ldmos_g3_same_bias_reclose import prepare


class TemplatesLdmosG3SameBiasRecloseTest(unittest.TestCase):
    @staticmethod
    def baseline() -> dict:
        return {
            "simulation_type": "sg_edge_flux_probe",
            "solver": {"mobility": {
                "model": "constant_field",
                "contact_electric_field_fallback": True,
            }},
            "contacts": [
                {"name": "gate", "bias": 0.0},
                {"name": "drain", "bias": 0.1},
            ],
            "sweep": {"bias_points": [0.0, 1.0]},
            "output_csv": "old.csv",
        }

    def test_prepare_preserves_controls_and_sets_gate(self) -> None:
        baseline = self.baseline()
        config = prepare(baseline, Path("state.csv"), 0.5, Path("out"))
        self.assertEqual(config["simulation_type"], "newton_solve_from_state")
        self.assertEqual(config["contacts"][0]["bias"], 0.5)
        self.assertEqual(config["contacts"][1]["bias"], 0.1)
        self.assertNotIn("sweep", config)
        self.assertNotIn("output_csv", config)
        self.assertIn("sweep", baseline)
        self.assertEqual(
            config["discretization"]["poisson_charge_volume_policy"],
            "global",
        )

    def test_prepare_enables_material_local_poisson_charge_volume(self) -> None:
        config = prepare(
            self.baseline(), Path("state.csv"), 1.0, Path("out"),
            "material_local")
        self.assertEqual(
            config["discretization"]["poisson_charge_volume_policy"],
            "material_local",
        )

    def test_prepare_fails_closed_on_unqualified_mobility(self) -> None:
        ialmob = self.baseline()
        ialmob["solver"]["mobility"]["ialmob"] = {"enabled": True}
        with self.assertRaisesRegex(ValueError, "IALMob disabled"):
            prepare(ialmob, Path("state.csv"), 0.5, Path("out"))

        missing_fallback = self.baseline()
        missing_fallback["solver"]["mobility"][
            "contact_electric_field_fallback"
        ] = False
        with self.assertRaisesRegex(ValueError, "contact HFS fallback"):
            prepare(missing_fallback, Path("state.csv"), 0.5, Path("out"))


if __name__ == "__main__":
    unittest.main()
