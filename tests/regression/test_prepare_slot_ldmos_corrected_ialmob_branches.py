import copy
import json
import tempfile
import unittest
from pathlib import Path

from scripts.prepare_slot_ldmos_corrected_ialmob_branches import (
    OUTPUT_ROOT,
    SEED_INNER_V,
    SEED_STATE,
    PreparationError,
    build_bootstrap_case,
    build_case,
    prepare,
    prepare_accelerated_extension,
    prepare_dense_low_voltage_extension,
    prepare_device_corrector_chunk,
    prepare_post_dense_extension,
    prepare_point_two_recovery,
    prepare_one_volt_extension,
    prepare_high_voltage_extension,
)


class CorrectedIalmobBranchPreparationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.base = {
            "output_csv": "outputs/stages/04_avalanche_activation_1v/iv.csv",
            "solver": {
                "verbose": True,
                "impact_ionization": {"source_jacobian": "frozen"},
                "mobility": {
                    "model": "masetti_field",
                    "doping_concentration_basis": "total_impurity",
                },
            },
            "sweep": {
                "start": 1.0,
                "stop": 1.0,
                "step": 1.0,
                "initial_state_file": "old.csv",
                "external_circuit": {"resistance_ohm_um": 1.0e12},
            },
        }

    def test_pair_has_one_controlled_physics_delta(self) -> None:
        off = build_case(copy.deepcopy(self.base), "ialmob_off", 16.0)
        on = build_case(copy.deepcopy(self.base), "ialmob_on", 16.0)
        self.assertEqual(off["solver"]["mobility"]["model"], "masetti_field")
        self.assertEqual(
            on["solver"]["mobility"]["model"], "masetti_field_lombardi"
        )
        self.assertEqual(
            on["solver"]["mobility"]["surface"]["surface_interface"],
            ["Silicon_1", "Oxide_1"],
        )
        self.assertNotIn("surface", off["solver"]["mobility"])
        self.assertEqual(
            off["solver"]["impact_ionization"]["source_jacobian"], "local_ad"
        )
        self.assertEqual(
            on["solver"]["impact_ionization"]["source_jacobian"], "local_ad"
        )
        self.assertFalse(off["sweep"]["external_circuit"]["enabled"])
        self.assertFalse(on["sweep"]["external_circuit"]["enabled"])
        self.assertEqual(off["solver"]["handoff"]["gummel_max_iter"], 0)
        self.assertEqual(on["solver"]["handoff"]["gummel_max_iter"], 0)

    def test_first_point_forces_reclosure_from_common_seed(self) -> None:
        bootstrap = build_bootstrap_case(
            copy.deepcopy(self.base), "ialmob_off"
        )
        self.assertEqual(bootstrap["sweep"]["initial_state_file"], SEED_STATE)
        self.assertEqual(bootstrap["sweep"]["bias_points"], [0.01])
        self.assertNotIn("impact_ionization", bootstrap["solver"])
        deck = build_case(copy.deepcopy(self.base), "ialmob_off", 1.0)
        self.assertEqual(deck["sweep"]["bias_points"][0], 0.011)
        self.assertNotEqual(deck["sweep"]["bias_points"][0], SEED_INNER_V)
        self.assertFalse(deck["sweep"]["boundary_control"]["resume"])

    def test_prepare_writes_pair_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / SEED_STATE).parent.mkdir(parents=True)
            (bundle / SEED_STATE).write_text("seed\n", encoding="utf-8")
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            manifest = prepare(bundle, 16.0)
            self.assertEqual(manifest["mobility_bootstrap_voltage_V"], 0.01)
            self.assertEqual(manifest["first_reclosed_voltage_V"], 0.011)
            self.assertEqual(manifest["bias_points_V"][-1], 16.0)
            for case in ("ialmob_off", "ialmob_on"):
                self.assertTrue(
                    (bundle / f"simulation_corrected_low_voltage_{case}.json").is_file()
                )
                self.assertTrue(
                    (bundle / f"simulation_corrected_low_voltage_{case}_bootstrap.json").is_file()
                )
                self.assertTrue((bundle / OUTPUT_ROOT / case / "states").is_dir())

    def test_rejects_missing_seed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            with self.assertRaisesRegex(PreparationError, "seed does not exist"):
                prepare(bundle, 1.0)

    def test_accelerated_extension_uses_secant_and_one_volt_cap(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            for case in ("ialmob_off", "ialmob_on"):
                seed = (
                    bundle / OUTPUT_ROOT / case / "states" /
                    "state_bias_0p050000.h5"
                )
                seed.parent.mkdir(parents=True)
                seed.write_text("seed\n", encoding="utf-8")
            manifest = prepare_accelerated_extension(bundle, 16.0)
            self.assertEqual(manifest["maximum_step_V"], 1.0)
            self.assertEqual(manifest["predictor"], "secant")
            for case in ("ialmob_off", "ialmob_on"):
                deck = json.loads(
                    (bundle / f"simulation_corrected_accelerated_{case}.json")
                    .read_text(encoding="utf-8")
                )
                self.assertEqual(deck["sweep"]["max_step"], 1.0)
                self.assertEqual(
                    deck["sweep"]["continuation"]["predictor"]["mode"],
                    "secant",
                )

    def test_corrector_chunk_restarts_same_bias_without_relaxed_forcing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            for case in ("ialmob_off", "ialmob_on"):
                source = (
                    bundle / "outputs/ialmob_ablation/corrected_accelerated_20260823" /
                    case / "rejected_states/attempt_1_bias_0p100000_final.csv"
                )
                source.parent.mkdir(parents=True)
                source.write_text("state\n", encoding="utf-8")
            manifest = prepare_device_corrector_chunk(bundle, 0.1, 1)
            self.assertFalse(manifest["forcing_relaxed"])
            self.assertEqual(manifest["maximum_newton_iterations"], 10)
            for item in manifest["cases"]:
                deck = json.loads((bundle / item["config"]).read_text())
                self.assertEqual(deck["sweep"]["bias_points"], [0.1])
                self.assertEqual(deck["solver"]["max_iter"], 10)
                self.assertEqual(
                    deck["solver"]["line_search_mode"], "block_filter"
                )
                self.assertEqual(
                    deck["solver"]["quasi_fermi_update_limit_mode"],
                    "uniform_trust_region",
                )
                self.assertEqual(
                    deck["solver"]["quasi_fermi_trust_region_growth_factor"], 2.0
                )
                self.assertEqual(
                    deck["solver"]["quasi_fermi_trust_region_max_multiplier"], 4.0
                )
                self.assertEqual(
                    deck["solver"]["quasi_fermi_trust_region_shrink_factor"], 0.5
                )
                self.assertEqual(
                    deck["solver"]["quasi_fermi_trust_region_min_multiplier"], 0.125
                )
                case = item["case"]
                rejected = (
                    bundle / "outputs/ialmob_ablation/device_manifold_corrector_20260823" /
                    "bias_0p100000" / case / "chunk_01/rejected_states"
                )
                self.assertTrue(rejected.is_dir())
            self.assertEqual(
                manifest["quasi_fermi_update_limit_mode"],
                "uniform_trust_region",
            )

    def test_dense_extension_uses_paired_seed_and_exact_10_mV_points(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            for case in ("ialmob_off", "ialmob_on"):
                seed = (
                    bundle / "outputs/ialmob_ablation/corrected_low_voltage_newton_20260823" /
                    case / "states/state_bias_0p050000.h5"
                )
                seed.parent.mkdir(parents=True)
                seed.write_text("state\n", encoding="utf-8")
            manifest = prepare_dense_low_voltage_extension(bundle, 0.1)
            self.assertEqual(manifest["bias_points_V"], [0.06, 0.07, 0.08, 0.09, 0.1])
            self.assertFalse(manifest["forcing_relaxed"])
            for item in manifest["cases"]:
                deck = json.loads((bundle / item["config"]).read_text())
                self.assertEqual(deck["sweep"]["max_step"], 0.01)
                self.assertEqual(deck["solver"]["max_iter"], 20)
                self.assertEqual(
                    deck["solver"]["quasi_fermi_update_limit_mode"],
                    "uniform_trust_region",
                )

    def test_post_dense_extension_starts_from_paired_100_mV_states(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            for case in ("ialmob_off", "ialmob_on"):
                seed = (
                    bundle / "outputs/ialmob_ablation/corrected_dense_low_voltage_20260823" /
                    case / "states/state_bias_0p100000.h5"
                )
                seed.parent.mkdir(parents=True)
                seed.write_text("state\n", encoding="utf-8")
            manifest = prepare_post_dense_extension(bundle, 0.5)
            self.assertEqual(
                manifest["bias_points_V"], [0.12, 0.15, 0.2, 0.3, 0.5]
            )
            self.assertFalse(manifest["forcing_relaxed"])
            for item in manifest["cases"]:
                deck = json.loads((bundle / item["config"]).read_text())
                self.assertEqual(deck["sweep"]["initial_step"], 0.02)
                self.assertEqual(deck["sweep"]["max_step"], 0.25)

    def test_point_two_recovery_disables_predictor_and_uses_10_mV_grid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            for case in ("ialmob_off", "ialmob_on"):
                seed = (
                    bundle / "outputs/ialmob_ablation/corrected_post_dense_20260823" /
                    case / "states/state_bias_0p150000.h5"
                )
                seed.parent.mkdir(parents=True)
                seed.write_text("state\n", encoding="utf-8")
            manifest = prepare_point_two_recovery(bundle)
            self.assertEqual(
                manifest["bias_points_V"], [0.16, 0.17, 0.18, 0.19, 0.2]
            )
            self.assertEqual(manifest["predictor"], "none")
            self.assertFalse(manifest["forcing_relaxed"])
            for item in manifest["cases"]:
                deck = json.loads((bundle / item["config"]).read_text())
                self.assertEqual(deck["sweep"]["max_step"], 0.01)
                self.assertEqual(
                    deck["sweep"]["continuation"]["predictor"]["mode"], "none"
                )
                self.assertEqual(
                    deck["solver"]["quasi_fermi_update_limit_mode"],
                    "uniform_trust_region",
                )

    def test_one_volt_extension_uses_paired_200_mV_anchors(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            for case in ("ialmob_off", "ialmob_on"):
                seed = (
                    bundle /
                    "outputs/ialmob_ablation/corrected_point_two_recovery_20260823" /
                    case / "states/state_bias_0p200000.h5"
                )
                seed.parent.mkdir(parents=True)
                seed.write_text("state\n", encoding="utf-8")
            manifest = prepare_one_volt_extension(bundle)
            self.assertEqual(
                manifest["bias_points_V"], [0.25, 0.3, 0.4, 0.5, 0.75, 1.0]
            )
            self.assertEqual(manifest["predictor"], "none")
            for item in manifest["cases"]:
                deck = json.loads((bundle / item["config"]).read_text())
                self.assertEqual(deck["sweep"]["initial_step"], 0.05)
                self.assertEqual(deck["sweep"]["max_step"], 0.25)
                self.assertEqual(
                    deck["sweep"]["continuation"]["predictor"]["mode"], "none"
                )

    def test_high_voltage_extension_uses_device_manifold_steps(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle / "simulation_04_avalanche_activation_1v.json").write_text(
                json.dumps(self.base), encoding="utf-8"
            )
            for case in ("ialmob_off", "ialmob_on"):
                seed = (
                    bundle /
                    "outputs/ialmob_ablation/corrected_one_volt_extension_20260823" /
                    case / "states/state_bias_1p000000.h5"
                )
                seed.parent.mkdir(parents=True)
                seed.write_text("state\n", encoding="utf-8")
            manifest = prepare_high_voltage_extension(bundle)
            self.assertEqual(manifest["seed_voltage_V"], 1.0)
            self.assertEqual(manifest["bias_points_V"][-1], 12.0)
            self.assertEqual(manifest["predictor"], "none")
            for item in manifest["cases"]:
                deck = json.loads((bundle / item["config"]).read_text())
                predictor = deck["sweep"]["continuation"]["predictor"]
                self.assertEqual(predictor["mode"], "none")
                self.assertEqual(deck["sweep"]["bias_points"][:4], [1.05, 1.1, 1.2, 1.3])
                self.assertEqual(deck["sweep"]["initial_step"], 0.05)
                self.assertEqual(deck["sweep"]["growth_factor"], 1.25)
                self.assertEqual(deck["sweep"]["max_step"], 0.25)


if __name__ == "__main__":
    unittest.main()
