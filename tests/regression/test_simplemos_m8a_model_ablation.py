from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "simplemos_m8a",
    ROOT / "scripts" / "run_simplemos_m8a_model_ablation.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
CONFIRMATION_SPEC = importlib.util.spec_from_file_location(
    "simplemos_m8a_confirmation",
    ROOT / "scripts" / "run_simplemos_m8a_confirmation.py")
assert CONFIRMATION_SPEC and CONFIRMATION_SPEC.loader
CONFIRMATION = importlib.util.module_from_spec(CONFIRMATION_SPEC)
CONFIRMATION_SPEC.loader.exec_module(CONFIRMATION)
CONTRACT_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8a_model_ablation_contract_v1.json"
)
CONFIRMATION_CONTRACT_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8a_confirmation_contract_v1.json"
)
EVIDENCE_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8a_model_ablation_evidence.json"
)


class SimpleMosM8AModelAblationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = MODULE.read_json(CONTRACT_PATH)
        self.variants = {
            item["id"]: item for item in self.contract["variants"]
        }

    def base_config(self) -> dict:
        return MODULE.m4.base_config(
            Path("mesh.json"), Path("doping.csv"), Path("materials.json"),
            {"model": "phumob_field_lombardi"})

    def test_contract_freezes_first_round_and_exact_bias_lattice(self) -> None:
        MODULE.validate_contract(self.contract)
        self.assertEqual("n23", self.contract["device"]["id"])
        self.assertEqual([0.05], self.contract["bias_matrix"][
            "drain_voltages_V"])
        self.assertEqual(
            [round(index * 0.05, 12) for index in range(51)],
            self.contract["bias_matrix"]["gate_lattice"]["values_V"])
        self.assertEqual("forbidden", self.contract["comparison"][
            "interpolation"])
        self.assertEqual(
            ["full", "no_srh", "plain_srh", "no_bgn", "no_enormal",
             "no_hfs", "phumob_only", "constant_mu"],
            list(self.variants))

    def test_sentaurus_decks_toggle_one_declared_model_group(self) -> None:
        full = MODULE.sentaurus_physics(self.variants["full"])
        self.assertIn("EffectiveIntrinsicDensity(OldSlotboom)", full)
        self.assertIn("Mobility(PhuMob HighFieldSaturation Enormal)", full)
        self.assertIn("Recombination(SRH(DopingDependence))", full)

        no_srh = MODULE.sentaurus_physics(self.variants["no_srh"])
        self.assertNotIn("Recombination", no_srh)
        self.assertIn("Mobility(PhuMob HighFieldSaturation Enormal)", no_srh)

        plain_srh = MODULE.sentaurus_physics(self.variants["plain_srh"])
        self.assertIn("Recombination(SRH)", plain_srh)
        self.assertNotIn("DopingDependence", plain_srh)

        no_bgn = MODULE.sentaurus_physics(self.variants["no_bgn"])
        self.assertNotIn("OldSlotboom", no_bgn)
        self.assertIn("Recombination(SRH(DopingDependence))", no_bgn)

        self.assertIn(
            "Mobility(PhuMob HighFieldSaturation)",
            MODULE.sentaurus_physics(self.variants["no_enormal"]))
        self.assertIn(
            "Mobility(PhuMob Enormal)",
            MODULE.sentaurus_physics(self.variants["no_hfs"]))
        self.assertIn(
            "Mobility(PhuMob)",
            MODULE.sentaurus_physics(self.variants["phumob_only"]))
        self.assertNotIn(
            "Mobility(",
            MODULE.sentaurus_physics(self.variants["constant_mu"]))

    def test_vela_mobility_mapping_matches_sentaurus_ablation(self) -> None:
        expected = {
            "full": "phumob_field_lombardi",
            "no_srh": "phumob_field_lombardi",
            "plain_srh": "phumob_field_lombardi",
            "no_bgn": "phumob_field_lombardi",
            "no_enormal": "phumob_field",
            "no_hfs": "phumob_lombardi",
            "phumob_only": "phumob",
            "constant_mu": "constant",
        }
        for variant_id, model in expected.items():
            config = self.base_config()
            MODULE.apply_variant_physics(
                config, self.variants[variant_id], self.contract)
            self.assertEqual(model, config["solver"]["mobility"]["model"])

    def test_vela_recombination_and_bgn_mappings_are_explicit(self) -> None:
        config = self.base_config()
        MODULE.apply_variant_physics(
            config, self.variants["no_srh"], self.contract)
        self.assertEqual(["none"], config["solver"]["recombination"])
        self.assertFalse(config["solver"]["srh_doping_dependence"]["enabled"])

        config = self.base_config()
        MODULE.apply_variant_physics(
            config, self.variants["plain_srh"], self.contract)
        self.assertEqual(["srh"], config["solver"]["recombination"])
        self.assertFalse(config["solver"]["srh_doping_dependence"]["enabled"])
        self.assertEqual(1.0e-5, config["solver"]["taun"])
        self.assertEqual(3.0e-6, config["solver"]["taup"])

        config = self.base_config()
        MODULE.apply_variant_physics(
            config, self.variants["no_bgn"], self.contract)
        self.assertEqual("none", config["solver"]["bandgap_narrowing"]["model"])

        config = self.base_config()
        MODULE.apply_variant_physics(
            config, self.variants["full"], self.contract)
        srh = config["solver"]["srh_doping_dependence"]
        self.assertTrue(srh["enabled"])
        self.assertEqual(1.0e-5, srh["electron"]["tau_max_s"])
        self.assertEqual(3.0e-6, srh["hole"]["tau_max_s"])
        self.assertEqual(1.0e16, srh["electron"]["reference_doping_m3"])

    def test_parallel_executors_reject_zero_workers(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one"):
            MODULE.execute_vela(
                self.contract, Path("unused"), Path("runner"), jobs=0)

    def test_confirmation_contract_freezes_hfs_controls(self) -> None:
        contract = CONFIRMATION.read_json(CONFIRMATION_CONTRACT_PATH)
        CONFIRMATION.validate_contract(contract)
        self.assertEqual(["n17", "n21"], [
            item["id"] for item in contract["devices"]])
        self.assertEqual(["full", "no_hfs", "phumob_only"], [
            item["id"] for item in contract["variants"]])
        self.assertEqual([0.05, 1.0], contract["bias_matrix"][
            "drain_voltages_V"])

    def test_confirmation_sentaurus_deck_uses_shared_tdr_and_no_hfs(self) -> None:
        contract = CONFIRMATION.read_json(CONFIRMATION_CONTRACT_PATH)
        variant = next(item for item in contract["variants"]
                       if item["id"] == "no_hfs")
        deck = MODULE.sentaurus_deck("fixture", 1.0, variant, 50).replace(
            'Grid="input_fps.tdr"', 'Grid="../input_fps.tdr"')
        self.assertIn('Grid="../input_fps.tdr"', deck)
        self.assertIn("Mobility(PhuMob Enormal)", deck)
        self.assertNotIn("HighFieldSaturation", deck)
        self.assertIn("CurrentPlot(Time=(Range=(0 1) Intervals=50))", deck)

    def test_frozen_evidence_records_complete_paired_diagnosis(self) -> None:
        evidence = MODULE.read_json(EVIDENCE_PATH)
        self.assertEqual("complete", evidence["status"])
        self.assertEqual(1020, evidence["execution"][
            "total_direct_bias_points"])
        self.assertEqual("forbidden", evidence["execution"]["interpolation"])
        self.assertTrue(evidence["findings"][
            "no_hfs_improves_all_confirmation_conditions"])
        self.assertTrue(all(value < 0.0 for value in evidence["findings"][
            "confirmation_no_hfs_changes_dex"]))
        self.assertEqual(7, len(evidence["artifacts"]["figures"]))

        for key, expected_count in (("first_round_report", 8),
                                    ("confirmation_report", 12)):
            artifact = evidence["artifacts"][key]
            report_path = ROOT / artifact["path"]
            self.assertEqual(artifact["sha256"], MODULE.sha256(report_path))
            report = MODULE.read_json(report_path)
            self.assertEqual(expected_count, len(report["cases"]))
            for case in report["cases"]:
                comparison = ROOT / case["comparison_csv"]
                self.assertEqual(case["comparison_csv_sha256"],
                                 MODULE.sha256(comparison))

        for figure in evidence["artifacts"]["figures"]:
            path = ROOT / figure["path"]
            self.assertTrue(path.is_file())
            self.assertEqual(figure["sha256"], MODULE.sha256(path))
        for relative, expected_hash in evidence[
                "implementation_sha256"].items():
            self.assertEqual(expected_hash, MODULE.sha256(ROOT / relative))


if __name__ == "__main__":
    unittest.main()
