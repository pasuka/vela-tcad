import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[2] / "scripts" / "prepare_slot_ldmos_parameter_parity_ladder.py"
SPEC = importlib.util.spec_from_file_location("prepare_parameter_parity_ladder", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PrepareParameterParityLadderTest(unittest.TestCase):
    def test_pair_is_strict_and_avalanche_free(self) -> None:
        base = {
            "output_csv": "old.csv",
            "solver": {
                "mobility": {"model": "masetti_field_lombardi", "surface": {}},
                "impact_ionization": {"model": "van_overstraeten"},
                "handoff": {},
            },
            "sweep": {"external_circuit": {"enabled": True}},
        }
        off = MODULE.build_case(base, "ialmob_off", 1.0)
        on = MODULE.build_case(base, "ialmob_on", 1.0)
        self.assertNotIn("impact_ionization", off["solver"])
        self.assertNotIn("external_circuit", off["sweep"])
        self.assertEqual(off["solver"]["mobility"]["model"], "masetti_field")
        self.assertEqual(on["solver"]["mobility"]["model"], "masetti_field_lombardi")
        self.assertEqual(
            on["solver"]["mobility"]["surface"]["surface_interface"],
            ["Silicon_1", "Oxide_1"],
        )
        self.assertEqual(
            MODULE._normalized(off, "ialmob_off"),
            MODULE._normalized(on, "ialmob_on"),
        )

    def test_prepare_uses_common_equilibrium_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "outputs/stages/00_equilibrium").mkdir(parents=True)
            (bundle / MODULE.EQUILIBRIUM_STATE).write_text("node_id,psi,phin,phip\n")
            base = {
                "output_csv": "old.csv",
                "solver": {"mobility": {"model": "masetti_field"}, "handoff": {}},
                "sweep": {},
            }
            (bundle / "simulation_01_unit_resistor_1v.json").write_text(
                json.dumps(base), encoding="utf-8"
            )
            manifest = MODULE.prepare(bundle, 0.2)
            self.assertEqual(manifest["common_seed_state"], MODULE.EQUILIBRIUM_STATE)
            self.assertEqual(manifest["bias_points_V"][-1], 0.2)
            for case in MODULE.CASES:
                deck = json.loads(
                    (bundle / f"simulation_parameter_parity_avalanche_off_{case}.json").read_text()
                )
                self.assertEqual(deck["sweep"]["initial_state_file"], MODULE.EQUILIBRIUM_STATE)
                self.assertEqual(deck["sweep"]["continuation"]["predictor"]["mode"], "none")

    def test_direct_restart_recloses_one_common_high_voltage_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "diagnostics").mkdir()
            (bundle / "diagnostics/state_60v.csv").write_text(
                "node_id,psi,phin,phip\n", encoding="utf-8")
            base = {
                "output_csv": "old.csv",
                "solver": {"mobility": {"model": "masetti_field"}, "handoff": {}},
                "sweep": {},
            }
            (bundle / "simulation_01_unit_resistor_1v.json").write_text(
                json.dumps(base), encoding="utf-8")
            manifest = MODULE.prepare_direct_restart(
                bundle, "diagnostics/state_60v.csv", 59.99)
            self.assertEqual(manifest["bias_voltage_V"], 59.99)
            for case in MODULE.CASES:
                deck = json.loads(
                    (bundle / f"simulation_parameter_parity_direct_restart_{case}.json")
                    .read_text(encoding="utf-8"))
                self.assertEqual(deck["sweep"]["bias_points"], [59.99])
                self.assertEqual(
                    deck["sweep"]["initial_state_file"],
                    "diagnostics/state_60v.csv")
                self.assertNotIn("impact_ionization", deck["solver"])
                self.assertEqual(deck["solver"]["method"], "gummel_newton")
                self.assertEqual(deck["solver"]["handoff"]["gummel_max_iter"], 100)

    def test_direct_poisson_preserves_imported_quasi_fermi_potentials(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "diagnostics").mkdir()
            (bundle / "diagnostics/state_60v.csv").write_text(
                "node_id,psi,phin,phip\n", encoding="utf-8")
            base = {
                "output_csv": "old.csv",
                "solver": {"mobility": {"model": "masetti_field"}, "handoff": {}},
                "sweep": {},
            }
            (bundle / "simulation_01_unit_resistor_1v.json").write_text(
                json.dumps(base), encoding="utf-8")
            manifest = MODULE.prepare_direct_poisson_reclosure(
                bundle, "diagnostics/state_60v.csv", 59.99)
            deck = json.loads((bundle / manifest["config"]).read_text())
            self.assertEqual(deck["solver"]["method"], "poisson_only")
            self.assertEqual(
                deck["solver"]["poisson_only_quasi_fermi_initialization"],
                "preserve_initial")
            self.assertNotIn("handoff", deck["solver"])


if __name__ == "__main__":
    unittest.main()
